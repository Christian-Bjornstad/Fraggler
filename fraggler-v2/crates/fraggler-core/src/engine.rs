use std::collections::BTreeMap;

use thiserror::Error;

use crate::contract::{
    EngineMessage, EnginePayload, ProgressEvent, RunKind, RunRequest, RunStatus, RunSummary,
    WarningEvent, WarningSeverity,
};

#[derive(Debug, Error)]
pub enum EngineError {
    #[error("run request is missing at least one input path")]
    MissingInputs,
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

    sink.emit(EngineMessage::new(
        request.correlation_id,
        EnginePayload::Warning(WarningEvent {
            severity: WarningSeverity::Warn,
            code: "engine_not_implemented".to_owned(),
            message: "fraggler-core is scaffolded, but the Rust engine port has not started yet."
                .to_owned(),
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
