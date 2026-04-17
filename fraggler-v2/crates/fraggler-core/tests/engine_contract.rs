use camino::Utf8PathBuf;
use fraggler_core::{
    AnalysisKind, ContractVersion, EngineMessage, EnginePayload, EventSink, InputSpec, OutputSpec,
    RunKind, RunOptions, RunRequest, RunStatus, run_request,
};
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
    let request = sample_request();
    let mut sink = CollectingSink::default();
    let result = run_request(&request, &mut sink);
    assert!(
        result.is_err(),
        "stub engine should still return not implemented"
    );
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
