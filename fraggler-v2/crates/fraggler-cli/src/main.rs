use std::fs;
use std::io::{self, BufRead, Write};

use anyhow::{Context, Result, bail};
use camino::Utf8PathBuf;
use clap::{Args, Parser, Subcommand, ValueEnum};
use fraggler_core::{
    AnalysisKind, ContractVersion, EngineMessage, InputSpec, OutputSpec, PrimitiveAnalysisResult,
    RunKind, RunOptions, RunRequest, analyze_fsa_primitives, run_request,
};
use serde::{Deserialize, Serialize};
use tracing_subscriber::EnvFilter;
use uuid::Uuid;

#[derive(Debug, Parser)]
#[command(
    name = "fraggler-cli",
    version,
    about = "HemaFrag Diagnostics Rust CLI"
)]
struct Cli {
    #[arg(long, global = true, default_value = "info")]
    log_filter: String,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Analyze(CommandArgs),
    Qc(CommandArgs),
    ValidateFlt3(CommandArgs),
    BuildReport(CommandArgs),
    ServePrimitives,
}

#[derive(Debug, Clone, Args)]
struct CommandArgs {
    #[arg(long)]
    json_request: Option<Utf8PathBuf>,
    #[arg(long, value_enum)]
    analysis: Option<AnalysisArg>,
    #[arg(long = "input")]
    inputs: Vec<Utf8PathBuf>,
    #[arg(long)]
    output_dir: Option<Utf8PathBuf>,
    #[arg(long)]
    report_dir: Option<Utf8PathBuf>,
    #[arg(long)]
    artifacts_dir: Option<Utf8PathBuf>,
    #[arg(long)]
    max_workers: Option<usize>,
    #[arg(long, default_value_t = true)]
    deterministic: bool,
    #[arg(long, default_value_t = true)]
    compact_json: bool,
    #[arg(long, default_value_t = true)]
    open_reports_in_browser: bool,
    #[arg(long, default_value_t = false)]
    shadow_reference_python: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum AnalysisArg {
    Clonality,
    Flt3,
    General,
}

impl From<AnalysisArg> for AnalysisKind {
    fn from(value: AnalysisArg) -> Self {
        match value {
            AnalysisArg::Clonality => Self::Clonality,
            AnalysisArg::Flt3 => Self::Flt3,
            AnalysisArg::General => Self::General,
        }
    }
}

#[derive(Debug, Deserialize)]
struct PrimitiveWorkerRequest {
    input: Utf8PathBuf,
    analysis: Option<AnalysisArg>,
}

#[derive(Debug, Serialize)]
struct PrimitiveWorkerResponse {
    ok: bool,
    result: Option<PrimitiveAnalysisResult>,
    error: Option<String>,
}

#[derive(Default)]
struct JsonLineSink;

impl fraggler_core::EventSink for JsonLineSink {
    fn emit(&mut self, message: EngineMessage) -> fraggler_core::engine::EngineResult<()> {
        let line = serde_json::to_string(&message)
            .map_err(|err| fraggler_core::EngineError::Sink(err.to_string()))?;
        println!("{line}");
        Ok(())
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::new(cli.log_filter))
        .with_target(false)
        .compact()
        .init();

    let (run_kind, args) = match cli.command {
        Commands::Analyze(args) => (RunKind::Analyze, args),
        Commands::Qc(args) => (RunKind::Qc, args),
        Commands::ValidateFlt3(args) => (RunKind::ValidateFlt3, args),
        Commands::BuildReport(args) => (RunKind::BuildReport, args),
        Commands::ServePrimitives => return serve_primitives(),
    };

    let request = load_request(run_kind, args)?;
    let mut sink = JsonLineSink;
    match run_request(&request, &mut sink) {
        Ok(summary) => {
            tracing::info!(status = ?summary.status, "run completed");
            Ok(())
        }
        Err(err) => bail!("{err}"),
    }
}

fn serve_primitives() -> Result<()> {
    let stdin = io::stdin();
    let mut stdout = io::stdout().lock();

    for line_result in stdin.lock().lines() {
        let line = line_result.context("failed to read worker request line")?;
        if line.trim().is_empty() {
            continue;
        }

        let response = match serde_json::from_str::<PrimitiveWorkerRequest>(&line) {
            Ok(request) => match analyze_fsa_primitives(
                &request.input,
                request.analysis.map(Into::into).as_ref(),
            ) {
                Ok(result) => PrimitiveWorkerResponse {
                    ok: true,
                    result: Some(result),
                    error: None,
                },
                Err(err) => PrimitiveWorkerResponse {
                    ok: false,
                    result: None,
                    error: Some(err.to_string()),
                },
            },
            Err(err) => PrimitiveWorkerResponse {
                ok: false,
                result: None,
                error: Some(format!("invalid worker request: {err}")),
            },
        };

        serde_json::to_writer(&mut stdout, &response)
            .context("failed to serialize worker response")?;
        stdout.write_all(b"\n").context("failed to write newline")?;
        stdout.flush().context("failed to flush worker response")?;
    }

    Ok(())
}

fn load_request(run_kind: RunKind, args: CommandArgs) -> Result<RunRequest> {
    if let Some(path) = args.json_request {
        let raw = fs::read_to_string(&path)
            .with_context(|| format!("failed to read JSON request from {path}"))?;
        let request: RunRequest = serde_json::from_str(&raw)
            .with_context(|| format!("failed to parse JSON request from {path}"))?;
        return Ok(request);
    }

    let output_dir = args
        .output_dir
        .context("--output-dir is required when --json-request is not used")?;

    Ok(RunRequest {
        contract_version: ContractVersion::default(),
        run_kind,
        analysis_kind: args.analysis.map(Into::into),
        correlation_id: Uuid::new_v4(),
        inputs: InputSpec {
            paths: args.inputs,
            manifest_path: None,
            report_source_path: None,
        },
        output: OutputSpec {
            root_dir: output_dir,
            report_dir: args.report_dir,
            artifacts_dir: args.artifacts_dir,
        },
        options: RunOptions {
            max_workers: args.max_workers,
            deterministic: args.deterministic,
            emit_compact_json: args.compact_json,
            open_reports_in_browser: args.open_reports_in_browser,
            shadow_reference_python: args.shadow_reference_python,
            ..RunOptions::default()
        },
    })
}
