use camino::Utf8PathBuf;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HtmlReportPlan {
    pub patient_id: String,
    pub output_path: Utf8PathBuf,
    pub title: String,
    pub includes_qc_section: bool,
    pub bundled_plotly: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReportPayload {
    pub patient_id: String,
    pub html_title: String,
    pub sections: Vec<ReportSection>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReportSection {
    pub id: String,
    pub heading: String,
    pub html_fragment: String,
}
