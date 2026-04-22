use std::path::{Path, PathBuf};
use std::thread;
use std::rc::Rc;

use anyhow::Result;
use camino::Utf8PathBuf;
use fraggler_core::{
    AnalysisKind, ContractVersion, EngineMessage, EnginePayload, EventSink, InputSpec, OutputSpec,
    RunKind, RunOptions, RunRequest, run_request,
};
use slint::{ComponentHandle, SharedString, Weak, Model, VecModel, ModelRc};
use uuid::Uuid;
use walkdir::WalkDir;

slint::include_modules!();

fn main() -> Result<()> {
    let app = MainWindow::new()?;
    let readme_path = workspace_readme_path();
    let default_output = default_output_path();

    // Setup input sources model
    let input_sources_model = Rc::new(VecModel::<SharedString>::default());
    app.set_input_sources(ModelRc::from(input_sources_model.clone()));

    app.set_output_path(SharedString::from(
        default_output.to_string_lossy().to_string(),
    ));
    app.set_last_artifact_path(SharedString::from(""));
    app.set_last_run_summary(SharedString::from("No runs yet."));
    app.set_active_hint(SharedString::from(analysis_hint("clonality")));

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
                app.set_active_hint(SharedString::from(analysis_hint(&analysis)));
            }
            append_log(&weak, &format!("[desktop] analysis set to {analysis}"));
            set_status(&weak, &format!("Selected analysis: {analysis}"));
        });
    }

    {
        let weak = app.as_weak();
        let model = input_sources_model.clone();
        app.on_add_files(move || {
            if let Some(paths) = rfd::FileDialog::new()
                .add_filter("FSA/ABI", &["fsa", "FSA", "abi", "ABI"])
                .pick_files()
            {
                for path in paths {
                    model.push(SharedString::from(path.to_string_lossy().to_string()));
                    append_log(
                        &weak,
                        &format!("[desktop] added file {}", path.display()),
                    );
                }
                set_status(&weak, "Files added to queue.");
            }
        });
    }

    {
        let weak = app.as_weak();
        let model = input_sources_model.clone();
        app.on_add_folder(move || {
            if let Some(path) = rfd::FileDialog::new().pick_folder() {
                model.push(SharedString::from(path.to_string_lossy().to_string()));
                append_log(
                    &weak,
                    &format!("[desktop] added folder {}", path.display()),
                );
                set_status(&weak, "Folder added to queue.");
            }
        });
    }

    {
        let model = input_sources_model.clone();
        app.on_remove_source(move |index| {
            model.remove(index as usize);
        });
    }

    {
        let weak = app.as_weak();
        app.on_browse_output(move || {
            if let Some(path) = rfd::FileDialog::new().pick_folder() {
                set_output_path(&weak, &path.to_string_lossy());
                append_log(
                    &weak,
                    &format!("[desktop] selected output {}", path.display()),
                );
                set_status(&weak, "Output directory selected.");
            }
        });
    }

    {
        let weak = app.as_weak();
        let sources_model = input_sources_model.clone();
        app.on_run_analysis(move || {
            let Some(app) = weak.upgrade() else {
                return;
            };

            let output_path = app.get_output_path().to_string();
            let analysis_kind = app.get_analysis_kind().to_string();

            let sources: Vec<String> = sources_model.iter().map(|s| s.to_string()).collect();

            if sources.is_empty() {
                set_status(&weak, "No input sources provided.");
                append_log(
                    &weak,
                    "[error] please add files or folders before running.",
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

            // Expand inputs
            let mut all_paths = Vec::new();
            for src in sources {
                let p = Path::new(&src);
                if p.is_dir() {
                    for entry in WalkDir::new(p).into_iter().filter_map(|e| e.ok()) {
                        if entry.file_type().is_file() {
                            let ext = entry.path().extension().and_then(|s| s.to_str()).unwrap_or("");
                            if ext.to_lowercase() == "fsa" {
                                all_paths.push(Utf8PathBuf::from(entry.path().to_string_lossy().to_string()));
                            }
                        }
                    }
                } else {
                    all_paths.push(Utf8PathBuf::from(src));
                }
            }

            if all_paths.is_empty() {
                set_status(&weak, "No .fsa files found in selected sources.");
                return;
            }

            let request = RunRequest {
                contract_version: ContractVersion::current(),
                run_kind: RunKind::Analyze,
                analysis_kind: Some(parse_analysis_kind(&analysis_kind)),
                correlation_id: Uuid::new_v4(),
                inputs: InputSpec {
                    paths: all_paths.clone(),
                    manifest_path: None,
                    report_source_path: None,
                },
                output: OutputSpec {
                    root_dir: Utf8PathBuf::from(output_path.clone()),
                    report_dir: None,
                    artifacts_dir: None,
                },
                options: RunOptions {
                    deterministic: true,
                    ..RunOptions::default()
                },
            };

            set_status(
                &weak,
                &format!("Running Rust analysis for {} file(s)", all_paths.len()),
            );
            set_last_artifact_path(&weak, "");
            set_last_run_summary(
                &weak,
                &format!("Processing {} files...", all_paths.len()),
            );
            append_log(
                &weak,
                &format!(
                    "[run] starting batch analyze | analysis={} | inputs={} | output={}",
                    analysis_kind, all_paths.len(), output_path
                ),
            );

            let weak_for_thread = weak.clone();
            thread::spawn(move || {
                let mut sink = DesktopEventSink::new(weak_for_thread.clone());
                match run_request(&request, &mut sink) {
                    Ok(summary) => {
                        let summary_path = Path::new(&output_path).join("analyze_summary.json");
                        if summary_path.exists() {
                            append_log(&weak_for_thread, "[report] triggering Python report builder bridge...");
                            match run_report_bridge(&summary_path, Path::new(&output_path)) {
                                Ok(_) => append_log(&weak_for_thread, "[report] Plotly reports generated and aggregated by DIT."),
                                Err(err) => append_log(&weak_for_thread, &format!("[report] bridge failed: {err}")),
                            }
                        }

                        set_status(
                            &weak_for_thread,
                            &format!("Rust analysis finished with status {:?}", summary.status),
                        );
                        set_last_run_summary(
                            &weak_for_thread,
                            &format!(
                                "Batch complete | {} | {} artifacts",
                                format!("{:?}", summary.status).to_lowercase(),
                                summary.artifact_manifest.len()
                            ),
                        );
                    }
                    Err(err) => {
                        append_log(&weak_for_thread, &format!("[error] {err}"));
                        set_status(&weak_for_thread, &format!("Run failed: {err}"));
                        set_last_run_summary(
                            &weak_for_thread,
                            &format!("Batch failed | {}", err),
                        );
                    }
                }
            });
        });
    }

    {
        let weak = app.as_weak();
        app.on_clear_log(move || {
            if let Some(app) = weak.upgrade() {
                app.set_log_text(SharedString::from(
                    "Rust + Slint desktop shell connected to fraggler-core.\n",
                ));
                set_status(&weak, "Log cleared.");
            }
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

fn set_last_run_summary(weak: &Weak<MainWindow>, summary: &str) {
    let weak = weak.clone();
    let summary = summary.to_owned();
    let _ = slint::invoke_from_event_loop(move || {
        if let Some(app) = weak.upgrade() {
            app.set_last_run_summary(SharedString::from(summary));
        }
    });
}

fn set_output_path(weak: &Weak<MainWindow>, path: &str) {
    let weak = weak.clone();
    let path = path.to_owned();
    let _ = slint::invoke_from_event_loop(move || {
        if let Some(app) = weak.upgrade() {
            app.set_output_path(SharedString::from(path));
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

fn analysis_hint(value: &str) -> &'static str {
    match value.to_ascii_lowercase().as_str() {
        "flt3" => "FLT3: use FLT3-ITD / D835 / NPM1-labelled .fsa files when available.",
        "general" => "General: useful for smoke tests and neutral signal exploration.",
        _ => "Clonality: use .fsa from IGK / KDE / TCRg-style workflows.",
    }
}

fn run_report_bridge(summary_path: &Path, output_dir: &Path) -> Result<()> {
    let root = workspace_root();
    let bridge_script = root.join("scripts").join("build_reports_v2.py");
    let python_exe = root.join("fraggler-mac310-venv").join("bin").join("python3");

    let mut cmd = std::process::Command::new(if python_exe.exists() {
        python_exe.as_os_str()
    } else {
        std::ffi::OsStr::new("python3")
    });

    cmd.arg(bridge_script)
       .arg(summary_path)
       .arg(output_dir);

    tracing::info!(?cmd, "running report bridge");
    let status = cmd.status()?;
    if !status.success() {
        anyhow::bail!("report bridge exited with error code {:?}", status.code());
    }
    Ok(())
}

fn workspace_root() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .map(|path| path.to_owned())
        .unwrap_or_else(|| PathBuf::from("."))
}

fn workspace_readme_path() -> PathBuf {
    workspace_root().join("README.md")
}

fn default_output_path() -> PathBuf {
    workspace_root()
        .join("validation_outputs")
        .join("fraggler_v2_desktop")
}

