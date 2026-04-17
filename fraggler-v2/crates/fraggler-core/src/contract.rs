use std::collections::BTreeMap;

use camino::Utf8PathBuf;
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ContractVersion {
    pub major: u16,
    pub minor: u16,
}

impl Default for ContractVersion {
    fn default() -> Self {
        Self { major: 1, minor: 0 }
    }
}

impl ContractVersion {
    pub fn current() -> Self {
        Self::default()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RunKind {
    Analyze,
    Qc,
    ValidateFlt3,
    BuildReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AnalysisKind {
    Clonality,
    Flt3,
    General,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct InputSpec {
    pub paths: Vec<Utf8PathBuf>,
    pub manifest_path: Option<Utf8PathBuf>,
    pub report_source_path: Option<Utf8PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct OutputSpec {
    pub root_dir: Utf8PathBuf,
    pub report_dir: Option<Utf8PathBuf>,
    pub artifacts_dir: Option<Utf8PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RunOptions {
    pub max_workers: Option<usize>,
    pub deterministic: bool,
    pub emit_compact_json: bool,
    pub open_reports_in_browser: bool,
    pub shadow_reference_python: bool,
    #[serde(default)]
    pub extra: BTreeMap<String, Value>,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            max_workers: None,
            deterministic: true,
            emit_compact_json: true,
            open_reports_in_browser: true,
            shadow_reference_python: false,
            extra: BTreeMap::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RunRequest {
    pub contract_version: ContractVersion,
    pub run_kind: RunKind,
    pub analysis_kind: Option<AnalysisKind>,
    pub correlation_id: Uuid,
    pub inputs: InputSpec,
    pub output: OutputSpec,
    pub options: RunOptions,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    Accepted,
    Running,
    Succeeded,
    Failed,
    NotImplemented,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProgressEvent {
    pub phase: String,
    pub file: Option<String>,
    pub files_done: Option<usize>,
    pub files_total: Option<usize>,
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WarningSeverity {
    Info,
    Warn,
    Error,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct WarningEvent {
    pub severity: WarningSeverity,
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReportArtifact {
    pub path: Utf8PathBuf,
    pub kind: String,
    pub size_bytes: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RunSummary {
    pub status: RunStatus,
    pub timings_ms: BTreeMap<String, u64>,
    pub artifact_manifest: Vec<ReportArtifact>,
    pub details: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "type", content = "data", rename_all = "snake_case")]
pub enum EnginePayload {
    RequestAccepted(RunRequest),
    Progress(ProgressEvent),
    Warning(WarningEvent),
    Artifact(ReportArtifact),
    Summary(RunSummary),
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct EngineMessage {
    pub timestamp: DateTime<Utc>,
    pub correlation_id: Uuid,
    pub payload: EnginePayload,
}

impl EngineMessage {
    pub fn new(correlation_id: Uuid, payload: EnginePayload) -> Self {
        Self {
            timestamp: Utc::now(),
            correlation_id,
            payload,
        }
    }
}
