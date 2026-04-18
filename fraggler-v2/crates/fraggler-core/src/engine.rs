use std::collections::BTreeMap;
use std::fs;

use camino::Utf8PathBuf;
use thiserror::Error;

use crate::contract::{
    EngineMessage, EnginePayload, ProgressEvent, ReportArtifact, RunKind, RunRequest, RunStatus,
    RunSummary, WarningEvent, WarningSeverity,
};
use crate::primitives::{PrimitiveAnalysisResult, analyze_fsa_primitives};

#[derive(Debug, Error)]
pub enum EngineError {
    #[error("run request is missing at least one input path")]
    MissingInputs,
    #[error("I/O error while trying to {context} at {path}: {source}")]
    Io {
        path: Utf8PathBuf,
        context: &'static str,
        #[source]
        source: std::io::Error,
    },
    #[error("invalid ABIF/FSA file at {path}: {message}")]
    InvalidAbif { path: Utf8PathBuf, message: String },
    #[error("signal math error: {0}")]
    SignalMath(String),
    #[error("primitive analysis failed: {message}")]
    PrimitiveAnalysis { message: String },
    #[error("run kind `{run_kind:?}` is not implemented yet in fraggler-core")]
    NotYetImplemented { run_kind: RunKind },
    #[error("event sink failed: {0}")]
    Sink(String),
}

pub type EngineResult<T> = Result<T, EngineError>;

pub trait EventSink {
    fn emit(&mut self, message: EngineMessage) -> EngineResult<()>;
}

#[derive(Default)]
pub struct NullSink;

impl EventSink for NullSink {
    fn emit(&mut self, _message: EngineMessage) -> EngineResult<()> {
        Ok(())
    }
}

pub fn run_request<S: EventSink>(request: &RunRequest, sink: &mut S) -> EngineResult<RunSummary> {
    if request.inputs.paths.is_empty() {
        return Err(EngineError::MissingInputs);
    }

    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::RequestAccepted(request.clone()),
    ))?;

    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Progress(ProgressEvent {
            phase: "bootstrap".to_owned(),
            file: None,
            files_done: Some(0),
            files_total: Some(request.inputs.paths.len()),
            note: Some("workspace skeleton accepted request".to_owned()),
        }),
    ))?;

    match request.run_kind {
        RunKind::Analyze => run_analyze_primitives(request, sink),
        _ => run_not_implemented(request, sink),
    }
}

fn run_analyze_primitives<S: EventSink>(
    request: &RunRequest,
    sink: &mut S,
) -> EngineResult<RunSummary> {
    let mut primitive_results = Vec::with_capacity(request.inputs.paths.len());
    for (index, input_path) in request.inputs.paths.iter().enumerate() {
        sink.emit(EngineMessage::new(
            request.correlation_id,
            EnginePayload::Progress(ProgressEvent {
                phase: "abif_parse".to_owned(),
                file: Some(input_path.to_string()),
                files_done: Some(index),
                files_total: Some(request.inputs.paths.len()),
                note: Some(
                    "phase 2 primitives: parsing ABIF and running native signal analysis"
                        .to_owned(),
                ),
            }),
        ))?;

        let primitive_result = analyze_fsa_primitives(input_path, request.analysis_kind.as_ref())?;
        primitive_results.push(primitive_result);
    }
    let artifact_manifest = persist_analyze_artifacts(request, &primitive_results, sink)?;
    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Warning(WarningEvent {
            severity: WarningSeverity::Warn,
            code: "partial_phase2_engine".to_owned(),
            message: "Native Rust ABIF parsing, baseline correction, ladder peak detection, and first-pass ladder candidate scoring are active; full ladder fitting is not ported yet."
                .to_owned(),
        }),
    ))?;

    let mut details = BTreeMap::from([
        (
            "run_kind".to_owned(),
            serde_json::to_value(&request.run_kind).unwrap_or_default(),
        ),
        (
            "analysis_kind".to_owned(),
            serde_json::to_value(&request.analysis_kind).unwrap_or_default(),
        ),
        (
            "primitive_stage".to_owned(),
            serde_json::Value::String("abif_baseline_peak_detection_ladder_candidates".to_owned()),
        ),
        (
            "primitive_results".to_owned(),
            serde_json::to_value(&primitive_results).unwrap_or_default(),
        ),
    ]);
    if let Some(first_result) = primitive_results.first() {
        details.insert(
            "primitive_result".to_owned(),
            serde_json::to_value(first_result).unwrap_or_default(),
        );
    }

    let summary = RunSummary {
        status: RunStatus::Succeeded,
        timings_ms: BTreeMap::from([("bootstrap".to_owned(), 0_u64)]),
        artifact_manifest,
        details,
    };

    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Summary(summary.clone()),
    ))?;
    Ok(summary)
}

fn persist_analyze_artifacts<S: EventSink>(
    request: &RunRequest,
    primitive_results: &[PrimitiveAnalysisResult],
    sink: &mut S,
) -> EngineResult<Vec<ReportArtifact>> {
    fs::create_dir_all(request.output.root_dir.as_std_path()).map_err(|source| {
        EngineError::Io {
            path: request.output.root_dir.clone(),
            context: "create output directory",
            source,
        }
    })?;

    let mut artifacts = Vec::new();

    let summary_path = request.output.root_dir.join("analyze_summary.json");
    let summary_json = serde_json::to_vec_pretty(primitive_results).map_err(|err| {
        EngineError::PrimitiveAnalysis {
            message: format!("failed to serialize analyze summary artifact: {err}"),
        }
    })?;
    fs::write(summary_path.as_std_path(), &summary_json).map_err(|source| EngineError::Io {
        path: summary_path.clone(),
        context: "write analyze summary artifact",
        source,
    })?;
    let summary_artifact = ReportArtifact {
        path: summary_path,
        kind: "analyze_summary_json".to_owned(),
        size_bytes: Some(summary_json.len() as u64),
    };
    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Artifact(summary_artifact.clone()),
    ))?;
    artifacts.push(summary_artifact);

    if let Some(first_result) = primitive_results.first() {
        let preview_path = request
            .output
            .root_dir
            .join("primitive_result_preview.json");
        let preview_json = serde_json::to_vec_pretty(first_result).map_err(|err| {
            EngineError::PrimitiveAnalysis {
                message: format!("failed to serialize primitive preview artifact: {err}"),
            }
        })?;
        fs::write(preview_path.as_std_path(), &preview_json).map_err(|source| EngineError::Io {
            path: preview_path.clone(),
            context: "write primitive preview artifact",
            source,
        })?;
        let preview_artifact = ReportArtifact {
            path: preview_path,
            kind: "primitive_result_preview_json".to_owned(),
            size_bytes: Some(preview_json.len() as u64),
        };
        sink.emit(EngineMessage::new(
            request.correlation_id,
            EnginePayload::Artifact(preview_artifact.clone()),
        ))?;
        artifacts.push(preview_artifact);
    }

    Ok(artifacts)
}

fn run_not_implemented<S: EventSink>(
    request: &RunRequest,
    sink: &mut S,
) -> EngineResult<RunSummary> {
    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Warning(WarningEvent {
            severity: WarningSeverity::Warn,
            code: "engine_not_implemented".to_owned(),
            message: "fraggler-core is scaffolded, but this run kind is not ported yet.".to_owned(),
        }),
    ))?;

    let summary = RunSummary {
        status: RunStatus::NotImplemented,
        timings_ms: BTreeMap::from([("bootstrap".to_owned(), 0_u64)]),
        artifact_manifest: Vec::new(),
        details: BTreeMap::from([
            (
                "run_kind".to_owned(),
                serde_json::to_value(&request.run_kind).unwrap_or_default(),
            ),
            (
                "analysis_kind".to_owned(),
                serde_json::to_value(&request.analysis_kind).unwrap_or_default(),
            ),
        ]),
    };

    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Summary(summary.clone()),
    ))?;

    Err(EngineError::NotYetImplemented {
        run_kind: request.run_kind.clone(),
    })
}
