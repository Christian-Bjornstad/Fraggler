use std::path::{Path, PathBuf};
use std::thread;

use anyhow::Result;
use camino::Utf8PathBuf;
use fraggler_core::{
    AnalysisKind, ContractVersion, EngineMessage, EnginePayload, EventSink, InputSpec, OutputSpec,
    RunKind, RunOptions, RunRequest, run_request,
};
use slint::{ComponentHandle, SharedString, Weak};
use uuid::Uuid;

slint::include_modules!();

fn main() -> Result<()> {
    let app = MainWindow::new()?;
    let readme_path = workspace_readme_path();
    let default_input = default_input_path();
    let default_output = default_output_path();

    app.set_input_path(SharedString::from(
        default_input.to_string_lossy().to_string(),
    ));
    app.set_output_path(SharedString::from(
        default_output.to_string_lossy().to_string(),
    ));
    app.set_last_artifact_path(SharedString::from(""));

    {
        let weak = app.as_weak();
        app.on_show_contract(move || {
            append_log(
                &weak,
                "[contract] v1 JSON contract active; desktop now runs fraggler-core directly.",
            );
            set_status(
                &weak,
                "Desktop shell is wired to the Rust engine. Ready for real analyze runs.",
            );
        });
    }

    {
        let weak = app.as_weak();
        app.on_select_analysis(move |analysis| {
            if let Some(app) = weak.upgrade() {
                app.set_analysis_kind(analysis.clone());
            }
            append_log(&weak, &format!("[desktop] analysis set to {analysis}"));
            set_status(&weak, &format!("Selected analysis: {analysis}"));
        });
    }

    {
        let weak = app.as_weak();
        app.on_run_analysis(move || {
            let Some(app) = weak.upgrade() else {
                return;
            };

            let input_path = app.get_input_path().to_string();
            let output_path = app.get_output_path().to_string();
            let analysis_kind = app.get_analysis_kind().to_string();

            if input_path.trim().is_empty() {
                set_status(&weak, "Input path is required.");
                append_log(
                    &weak,
                    "[error] please provide an input .fsa path before running.",
                );
                return;
            }
            if output_path.trim().is_empty() {
                set_status(&weak, "Output path is required.");
                append_log(
                    &weak,
                    "[error] please provide an output directory before running.",
                );
                return;
            }

            let request = RunRequest {
                contract_version: ContractVersion::current(),
                run_kind: RunKind::Analyze,
                analysis_kind: Some(parse_analysis_kind(&analysis_kind)),
                correlation_id: Uuid::new_v4(),
                inputs: InputSpec {
                    paths: vec![Utf8PathBuf::from(input_path.clone())],
                    manifest_path: None,
                    report_source_path: None,
                },
                output: OutputSpec {
                    root_dir: Utf8PathBuf::from(output_path.clone()),
                    report_dir: None,
                    artifacts_dir: None,
                },
                options: RunOptions::default(),
            };

            set_status(
                &weak,
                &format!("Running Rust analysis for {}", file_name(&input_path)),
            );
            set_last_artifact_path(&weak, "");
            append_log(
                &weak,
                &format!(
                    "[run] starting analyze | analysis={} | input={} | output={}",
                    analysis_kind, input_path, output_path
                ),
            );

            let weak_for_thread = weak.clone();
            thread::spawn(move || {
                let mut sink = DesktopEventSink::new(weak_for_thread.clone());
                match run_request(&request, &mut sink) {
                    Ok(summary) => {
                        let finished = format!(
                            "[done] status={:?} artifacts={} detail_keys={}",
                            summary.status,
                            summary.artifact_manifest.len(),
                            summary.details.len()
                        );
                        append_log(&weak_for_thread, &finished);
                        set_status(
                            &weak_for_thread,
                            &format!("Rust analysis finished with status {:?}", summary.status),
                        );
                    }
                    Err(err) => {
                        append_log(&weak_for_thread, &format!("[error] {err}"));
                        set_status(&weak_for_thread, &format!("Run failed: {err}"));
                    }
                }
            });
        });
    }

    {
        let weak = app.as_weak();
        app.on_open_output(move || {
            if let Some(app) = weak.upgrade() {
                let path = app.get_output_path().to_string();
                if let Err(err) = open::that(&path) {
                    append_log(&weak, &format!("[error] failed to open output path: {err}"));
                    set_status(&weak, &format!("Failed to open output path: {err}"));
                }
            }
        });
    }

    {
        let weak = app.as_weak();
        app.on_open_last_artifact(move || {
            if let Some(app) = weak.upgrade() {
                let path = app.get_last_artifact_path().to_string();
                if path.trim().is_empty() {
                    set_status(&weak, "No artifact has been produced yet.");
                    return;
                }
                if let Err(err) = open::that(&path) {
                    append_log(&weak, &format!("[error] failed to open artifact: {err}"));
                    set_status(&weak, &format!("Failed to open artifact: {err}"));
                }
            }
        });
    }

    app.on_open_docs(move || {
        if let Err(err) = open::that(&readme_path) {
            eprintln!("failed to open docs: {err}");
        }
    });

    app.run()?;
    Ok(())
}

struct DesktopEventSink {
    weak: Weak<MainWindow>,
}

impl DesktopEventSink {
    fn new(weak: Weak<MainWindow>) -> Self {
        Self { weak }
    }
}

impl EventSink for DesktopEventSink {
    fn emit(&mut self, message: EngineMessage) -> fraggler_core::engine::EngineResult<()> {
        let (status, line) = summarize_engine_message(&message);
        if let Some(status) = status {
            set_status(&self.weak, &status);
        }
        if let EnginePayload::Artifact(artifact) = &message.payload {
            set_last_artifact_path(&self.weak, artifact.path.as_str());
        }
        append_log(&self.weak, &line);
        Ok(())
    }
}

fn summarize_engine_message(message: &EngineMessage) -> (Option<String>, String) {
    match &message.payload {
        EnginePayload::RequestAccepted(request) => (
            Some(format!(
                "Accepted {} input(s) for {:?}",
                request.inputs.paths.len(),
                request.analysis_kind
            )),
            format!(
                "[accepted] run_kind={:?} analysis={:?} inputs={}",
                request.run_kind,
                request.analysis_kind,
                request.inputs.paths.len()
            ),
        ),
        EnginePayload::Progress(progress) => (
            Some(
                progress
                    .note
                    .clone()
                    .unwrap_or_else(|| format!("Phase {}", progress.phase)),
            ),
            format!(
                "[progress] phase={} file={} done={:?}/{:?} note={}",
                progress.phase,
                progress.file.clone().unwrap_or_else(|| "-".to_owned()),
                progress.files_done,
                progress.files_total,
                progress.note.clone().unwrap_or_else(|| "-".to_owned())
            ),
        ),
        EnginePayload::Warning(warning) => (
            Some(warning.message.clone()),
            format!(
                "[warning] severity={:?} code={} message={}",
                warning.severity, warning.code, warning.message
            ),
        ),
        EnginePayload::Artifact(artifact) => (
            Some(format!("Artifact emitted: {}", artifact.path)),
            format!(
                "[artifact] kind={} path={} size={:?}",
                artifact.kind, artifact.path, artifact.size_bytes
            ),
        ),
        EnginePayload::Summary(summary) => (
            Some(format!("Run summary: {:?}", summary.status)),
            format!(
                "[summary] status={:?} timings={:?} artifacts={} details={}",
                summary.status,
                summary.timings_ms,
                summary.artifact_manifest.len(),
                summary.details.len()
            ),
        ),
    }
}

fn append_log(weak: &Weak<MainWindow>, line: &str) {
    let weak = weak.clone();
    let line = line.to_owned();
    let _ = slint::invoke_from_event_loop(move || {
        if let Some(app) = weak.upgrade() {
            let mut log = app.get_log_text().to_string();
            if !log.ends_with('\n') {
                log.push('\n');
            }
            log.push_str(&line);
            log.push('\n');
            app.set_log_text(SharedString::from(log));
        }
    });
}

fn set_status(weak: &Weak<MainWindow>, status: &str) {
    let weak = weak.clone();
    let status = status.to_owned();
    let _ = slint::invoke_from_event_loop(move || {
        if let Some(app) = weak.upgrade() {
            app.set_status_text(SharedString::from(status));
        }
    });
}

fn set_last_artifact_path(weak: &Weak<MainWindow>, path: &str) {
    let weak = weak.clone();
    let path = path.to_owned();
    let _ = slint::invoke_from_event_loop(move || {
        if let Some(app) = weak.upgrade() {
            app.set_last_artifact_path(SharedString::from(path));
        }
    });
}

fn parse_analysis_kind(value: &str) -> AnalysisKind {
    match value.to_ascii_lowercase().as_str() {
        "flt3" => AnalysisKind::Flt3,
        "general" => AnalysisKind::General,
        _ => AnalysisKind::Clonality,
    }
}

fn file_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path)
        .to_owned()
}

fn workspace_readme_path() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .map(|path| path.join("README.md"))
        .unwrap_or_else(|| PathBuf::from("README.md"))
}

fn default_output_path() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .map(|path| path.join("validation_outputs").join("fraggler_v2_desktop"))
        .unwrap_or_else(|| PathBuf::from("validation_outputs/fraggler_v2_desktop"))
}

fn default_input_path() -> PathBuf {
    PathBuf::from(
        "/Volumes/T7 Shield/DATA/2026/2026_03_27_TCRg_IGK_KDE_CFB_H9H1DI2F_2026-03-27_0652/26OUM04817_IGK_270326_B05_H9H1DI2F.fsa",
    )
}
