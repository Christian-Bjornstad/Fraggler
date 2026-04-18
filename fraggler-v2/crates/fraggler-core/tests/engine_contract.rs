use camino::Utf8PathBuf;
use fraggler_core::{
    AnalysisKind, ContractVersion, EngineMessage, EnginePayload, EventSink, InputSpec, OutputSpec,
    RunKind, RunOptions, RunRequest, RunStatus, run_request,
};
use std::fs;
use std::path::Path;
use uuid::Uuid;

#[derive(Default)]
struct CollectingSink {
    messages: Vec<EngineMessage>,
}

impl EventSink for CollectingSink {
    fn emit(&mut self, message: EngineMessage) -> fraggler_core::engine::EngineResult<()> {
        self.messages.push(message);
        Ok(())
    }
}

fn sample_request() -> RunRequest {
    RunRequest {
        contract_version: ContractVersion::current(),
        run_kind: RunKind::Analyze,
        analysis_kind: Some(AnalysisKind::Clonality),
        correlation_id: Uuid::new_v4(),
        inputs: InputSpec {
            paths: vec![Utf8PathBuf::from("/tmp/sample.fsa")],
            manifest_path: None,
            report_source_path: None,
        },
        output: OutputSpec {
            root_dir: Utf8PathBuf::from("/tmp/out"),
            report_dir: Some(Utf8PathBuf::from("/tmp/out/REPORTS")),
            artifacts_dir: Some(Utf8PathBuf::from("/tmp/out/artifacts")),
        },
        options: RunOptions::default(),
    }
}

#[test]
fn contract_round_trip_serializes_cleanly() {
    let request = sample_request();
    let json = serde_json::to_string(&request).expect("serialize request");
    let decoded: RunRequest = serde_json::from_str(&json).expect("deserialize request");
    assert_eq!(decoded.run_kind, RunKind::Analyze);
    assert_eq!(decoded.analysis_kind, Some(AnalysisKind::Clonality));
    assert_eq!(decoded.inputs.paths.len(), 1);
    assert!(decoded.options.deterministic);
    assert!(decoded.options.emit_compact_json);
}

#[test]
fn run_request_emits_expected_stub_lifecycle() {
    let mut request = sample_request();
    request.run_kind = RunKind::Qc;
    let mut sink = CollectingSink::default();
    let result = run_request(&request, &mut sink);
    assert!(result.is_err(), "qc should still return not implemented");
    assert_eq!(sink.messages.len(), 4);

    assert!(matches!(
        sink.messages[0].payload,
        EnginePayload::RequestAccepted(_)
    ));
    assert!(matches!(
        sink.messages[1].payload,
        EnginePayload::Progress(_)
    ));
    assert!(matches!(
        sink.messages[2].payload,
        EnginePayload::Warning(_)
    ));

    match &sink.messages[3].payload {
        EnginePayload::Summary(summary) => {
            assert_eq!(summary.status, RunStatus::NotImplemented);
            assert!(summary.artifact_manifest.is_empty());
        }
        other => panic!("expected summary payload, got {other:?}"),
    }
}

#[test]
fn run_request_rejects_missing_inputs_without_emitting_events() {
    let mut request = sample_request();
    request.inputs.paths.clear();
    let mut sink = CollectingSink::default();
    let result = run_request(&request, &mut sink);
    assert!(result.is_err());
    assert!(sink.messages.is_empty());
}

#[test]
fn analyze_request_emits_bootstrap_before_io_failure() {
    let request = sample_request();
    let mut sink = CollectingSink::default();
    let result = run_request(&request, &mut sink);
    assert!(result.is_err());
    assert_eq!(sink.messages.len(), 3);

    assert!(matches!(
        sink.messages[0].payload,
        EnginePayload::RequestAccepted(_)
    ));
    assert!(matches!(
        sink.messages[1].payload,
        EnginePayload::Progress(_)
    ));
    assert!(matches!(
        sink.messages[2].payload,
        EnginePayload::Progress(_)
    ));
}

#[test]
fn analyze_request_writes_artifact_files_when_real_input_exists() {
    let real_input = Path::new(
        "/Volumes/T7 Shield/DATA/2026/2026_03_27_TCRg_IGK_KDE_CFB_H9H1DI2F_2026-03-27_0652/26OUM04817_IGK_270326_B05_H9H1DI2F.fsa",
    );
    if !real_input.exists() {
        return;
    }

    let temp_root =
        std::env::temp_dir().join(format!("fraggler_v2_artifact_test_{}", Uuid::new_v4()));
    let mut request = sample_request();
    request.inputs.paths = vec![Utf8PathBuf::from(real_input.to_string_lossy().to_string())];
    request.output.root_dir = Utf8PathBuf::from(temp_root.to_string_lossy().to_string());

    let mut sink = CollectingSink::default();
    let summary =
        run_request(&request, &mut sink).expect("analyze request should succeed for real input");
    assert_eq!(summary.status, RunStatus::Succeeded);
    assert!(!summary.artifact_manifest.is_empty());

    let artifact_paths = summary
        .artifact_manifest
        .iter()
        .map(|artifact| artifact.path.as_str().to_owned())
        .collect::<Vec<_>>();
    assert!(
        artifact_paths
            .iter()
            .any(|path| path.ends_with("analyze_summary.json"))
    );
    assert!(
        artifact_paths
            .iter()
            .any(|path| path.ends_with("primitive_result_preview.json"))
    );
    for path in artifact_paths {
        assert!(
            Path::new(&path).exists(),
            "expected artifact file {path} to exist"
        );
    }

    let _ = fs::remove_dir_all(temp_root);
}
