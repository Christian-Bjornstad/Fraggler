use std::collections::BTreeMap;

use camino::Utf8Path;
use serde::{Deserialize, Serialize};

use crate::abif::AbifRecord;
use crate::contract::AnalysisKind;
use crate::engine::EngineError;
use crate::ladders::LadderKind;
use crate::signal::{
    Peak, baseline_correct_guarded_nonnegative, baseline_correct_quantile_nonnegative, find_peaks,
};

const MAX_CANDIDATE_COMBINATIONS: usize = 2_000_000;
const BEAM_SEARCH_TRIGGER_COMBINATIONS: usize = 50_000;
const BEAM_SEARCH_WIDTH: usize = 192;
const BEAM_SEARCH_FINAL_CAP: usize = 4096;
const LADDER_MAX_GAP_EXPANSIONS: usize = 30;
const LADDER_GAP_EXPANSION_STEP: usize = 10;
const MAX_REFINEMENT_STEPS: usize = 3;
const MAX_REFINEMENT_OPTIONS_PER_STEP: usize = 5;
const MIN_REFINEMENT_TRIGGER_BP: f64 = 0.75;
const MAX_REFINEMENT_RADIUS_SCANS: f64 = 120.0;
const LADDER_DOMAIN_TIME_WEIGHT: f64 = 3.00;
const LADDER_DOMAIN_GAP_WEIGHT: f64 = 0.15;
const LADDER_DOMAIN_INTENSITY_WEIGHT: f64 = 0.10;
const LADDER_DOMAIN_FIRST_ANCHOR_WEIGHT: f64 = 1.50;
const ROX_HARD_TIME_MIN: f64 = 1300.0;
const ROX_HARD_TIME_MAX: f64 = 4300.0;
const ROX_PREFERRED_TIME_MIN: f64 = 1500.0;
const ROX_PREFERRED_TIME_MAX: f64 = 4000.0;
const ROX_MAX_FIRST_ANCHOR: f64 = 1900.0;
const GS500ROX_HARD_TIME_MIN: f64 = 1380.0;
const GS500ROX_HARD_TIME_MAX: f64 = 4550.0;
const GS500ROX_PREFERRED_TIME_MIN: f64 = 1490.0;
const GS500ROX_PREFERRED_TIME_MAX: f64 = 4425.0;
const GS500ROX_MAX_FIRST_ANCHOR: f64 = 1700.0;
const LIZ_HARD_TIME_MIN: f64 = 1150.0;
const LIZ_HARD_TIME_MAX: f64 = 4300.0;
const LIZ_PREFERRED_TIME_MIN: f64 = 1250.0;
const LIZ_PREFERRED_TIME_MAX: f64 = 4100.0;
const LIZ_MAX_FIRST_ANCHOR: f64 = 1700.0;
const SAMPLE_ASSAY_GROUP_DISTANCE_BP: f64 = 12.0;
const SAMPLE_ASSAY_MIN_RATIO: f64 = 0.40;
const CLONAL_MAX_LABELLED_PEAKS: usize = 3;
const CLONAL_DOMINANCE_RATIO: f64 = 1.7;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LadderFitPreview {
    pub max_allowed_peak_gap: usize,
    pub gap_expansions: usize,
    pub estimated_combination_count: usize,
    pub candidate_generation_capped: bool,
    pub evaluated_combination_count: usize,
    pub best_scan_indices: Vec<usize>,
    pub best_curvature_score: Option<f64>,
    pub best_quadratic_r2: Option<f64>,
    pub sizing_model: Option<SizingModelPreview>,
    pub refinement: Option<RefinementPreview>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SizingModelPreview {
    pub strategy: String,
    pub degree: usize,
    pub coefficients: Vec<f64>,
    pub predicted_ladder_basepairs: Vec<f64>,
    pub qc_metrics: LadderQcMetrics,
    pub sample_mapping: Option<SampleMappingPreview>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RefinementPreview {
    pub changed_step_indices: Vec<usize>,
    pub original_scan_indices: Vec<usize>,
    pub refined_scan_indices: Vec<usize>,
    pub refined_curvature_score: f64,
    pub refined_quadratic_r2: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SampleMappingPreview {
    pub points_retained: usize,
    pub min_basepair: f64,
    pub max_basepair: f64,
    pub monotonic_unique: bool,
    pub preview: Vec<SampleMappingPoint>,
    pub sample_peak_preview: Vec<SamplePeakPreview>,
    pub assay_group_preview: Vec<SamplePeakGroupPreview>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SampleMappingPoint {
    pub time: usize,
    pub intensity: f64,
    pub basepair: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SamplePeakPreview {
    pub time: usize,
    pub intensity: f64,
    pub basepair: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SamplePeakGroupPreview {
    pub group_id: usize,
    pub start_basepair: f64,
    pub end_basepair: f64,
    pub cluster_width_bp: f64,
    pub max_intensity: f64,
    pub dominant_peak_basepair: f64,
    pub dominant_peak_intensity: f64,
    pub dominant_ratio_vs_second: Option<f64>,
    pub kept_peak_count: usize,
    pub clonal_candidate: bool,
    pub peaks: Vec<SamplePeakPreview>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LadderQcMetrics {
    pub r2: f64,
    pub mean_abs_error_bp: f64,
    pub max_abs_error_bp: f64,
    pub monotonic_on_ladder: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PrimitiveAnalysisResult {
    pub file_name: String,
    pub scan_count: usize,
    pub data_channels: Vec<String>,
    pub dye_names: BTreeMap<String, String>,
    pub ladder: String,
    pub sample_channel_guess: String,
    pub size_standard_channel_guess: String,
    pub ladder_peak_count: usize,
    pub ladder_peak_preview: Vec<Peak>,
    pub ladder_fit_preview: Option<LadderFitPreview>,
    pub clonality_preview: Option<ClonalityPreview>,
    pub flt3_preview: Option<Flt3Preview>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ClonalityPreview {
    pub sample_channel: String,
    pub ranked_assays: Vec<ClonalityAssayMatch>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ClonalityAssayMatch {
    pub assay_name: String,
    pub matched_by_filename: bool,
    pub compatible_channel: bool,
    pub score: f64,
    pub clonal_group_count: usize,
    pub best_dominant_ratio: Option<f64>,
    pub matched_groups: Vec<ClonalityGroupMatch>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ClonalityGroupMatch {
    pub group_id: usize,
    pub overlap_start_bp: f64,
    pub overlap_end_bp: f64,
    pub peak_count: usize,
    pub cluster_width_bp: f64,
    pub dominant_peak_basepair: f64,
    pub dominant_peak_intensity: f64,
    pub dominant_ratio_vs_second: Option<f64>,
    pub clonal_candidate: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Flt3Preview {
    pub assay_name: String,
    pub matched_by_filename: bool,
    pub compatible_channel: bool,
    pub assay_bp_min: f64,
    pub assay_bp_max: f64,
    pub wt_peak: Option<SamplePeakPreview>,
    pub mutant_peaks: Vec<SamplePeakPreview>,
    pub strongest_mutant_ratio: Option<f64>,
    pub positive_call: bool,
}

#[derive(Debug, Clone, Copy)]
struct Flt3AssayDef {
    name: &'static str,
    channels: &'static [&'static str],
    bp_min: f64,
    bp_max: f64,
    wt_bp: f64,
    wt_tolerance_bp: f64,
    mutant_bp_min: f64,
    mutant_bp_max: f64,
    positive_ratio: f64,
    aliases: &'static [&'static str],
}

pub fn analyze_fsa_primitives(
    path: &Utf8Path,
    analysis_kind: Option<&AnalysisKind>,
) -> Result<PrimitiveAnalysisResult, EngineError> {
    let record = AbifRecord::from_path(path)?;
    let data_channels = record.data_channels();
    let size_standard_channel =
        select_size_standard_channel(&record, analysis_kind).ok_or_else(|| {
            EngineError::PrimitiveAnalysis {
                message: "no usable size-standard channel was found in the ABIF record".to_owned(),
            }
        })?;
    let file_name = record.file_name.clone();
    let sample_channel = preferred_sample_channel(
        &data_channels,
        &size_standard_channel,
        &file_name,
        analysis_kind,
    )
    .unwrap_or_else(|| {
        data_channels
            .iter()
            .find(|channel| channel.as_str() != size_standard_channel)
            .cloned()
            .unwrap_or_else(|| size_standard_channel.clone())
    });

    let ladder = suggested_ladder_kind(&record, &size_standard_channel, analysis_kind);
    let raw_min_height = match ladder {
        LadderKind::Liz500250 => 150.0,
        // ROX traces in clonality are prone to dense late noise clusters.
        // A stricter floor keeps true ladder peaks while suppressing blob tails.
        LadderKind::Rox400Hd => 180.0,
        LadderKind::Gs500Rox => 120.0,
    };
    // Baseline correction compresses peak amplitudes, so use a softer floor
    // on the corrected trace while still keeping obvious baseline chatter out.
    let corrected_min_height = (raw_min_height * 0.55_f64).max(50.0_f64);
    let min_distance = match ladder {
        LadderKind::Liz500250 => 15,
        LadderKind::Rox400Hd => 14,
        LadderKind::Gs500Rox => 15,
    };

    let size_standard_trace = record
        .channel_values(&size_standard_channel)
        .ok_or_else(|| EngineError::PrimitiveAnalysis {
            message: format!(
                "size-standard channel {size_standard_channel} is missing numeric data"
            ),
        })?;
    let sample_trace =
        record
            .channel_values(&sample_channel)
            .ok_or_else(|| EngineError::PrimitiveAnalysis {
                message: format!("sample channel {sample_channel} is missing numeric data"),
            })?;
    let corrected =
        baseline_correct_guarded_nonnegative(&size_standard_trace, 0.99, 100.0, 1000, 200, 0.10)?;
    let quantile_corrected =
        baseline_correct_quantile_nonnegative(&size_standard_trace, 200, 0.10);
    let ladder_peaks = select_ladder_peaks(
        &size_standard_trace,
        &corrected,
        &quantile_corrected,
        raw_min_height,
        corrected_min_height,
        min_distance,
        ladder.expected_peak_count() + 15,
        ladder.expected_peak_count(),
    );
    let dye_names = dye_names(&record);
    let ladder_fit_preview = build_ladder_fit_preview(&ladder_peaks, &sample_trace, ladder);
    let clonality_preview = if matches!(analysis_kind, Some(AnalysisKind::Clonality)) {
        ladder_fit_preview.as_ref().and_then(|preview| {
            preview
                .sizing_model
                .as_ref()
                .and_then(|model| model.sample_mapping.as_ref())
                .map(|mapping| {
                    build_clonality_preview(
                        &file_name,
                        &sample_channel,
                        &mapping.assay_group_preview,
                    )
                })
        })
    } else {
        None
    };
    let flt3_preview = if matches!(analysis_kind, Some(AnalysisKind::Flt3)) {
        ladder_fit_preview.as_ref().and_then(|preview| {
            preview
                .sizing_model
                .as_ref()
                .and_then(|model| model.sample_mapping.as_ref())
                .map(|mapping| {
                    build_flt3_preview(&file_name, &sample_channel, &mapping.sample_peak_preview)
                })
        })
    } else {
        None
    };

    Ok(PrimitiveAnalysisResult {
        file_name,
        scan_count: size_standard_trace.len(),
        data_channels,
        dye_names,
        ladder: ladder.display_name().to_owned(),
        sample_channel_guess: sample_channel,
        size_standard_channel_guess: size_standard_channel,
        ladder_peak_count: ladder_peaks.len(),
        ladder_peak_preview: ladder_peaks.into_iter().take(8).collect(),
        ladder_fit_preview,
        clonality_preview,
        flt3_preview,
    })
}

fn dye_names(record: &AbifRecord) -> BTreeMap<String, String> {
    let mut dyes = BTreeMap::new();
    for index in 1..=8 {
        let key = format!("DyeN{index}");
        if let Some(value) = record.string_value(&key) {
            dyes.insert(key, value.to_owned());
        }
    }
    dyes
}

fn select_size_standard_channel(
    record: &AbifRecord,
    analysis_kind: Option<&AnalysisKind>,
) -> Option<String> {
    match analysis_kind {
        Some(AnalysisKind::Clonality) => {
            if record.tags.contains_key("DATA105") {
                return Some("DATA105".to_owned());
            }
            if record.tags.contains_key("DATA4") {
                return Some("DATA4".to_owned());
            }
        }
        Some(AnalysisKind::Flt3) => {
            if record.tags.contains_key("DATA4") {
                return Some("DATA4".to_owned());
            }
            if record.tags.contains_key("DATA105") {
                return Some("DATA105".to_owned());
            }
        }
        _ => {}
    }

    for preferred in ["DATA105", "DATA4"] {
        if record.tags.contains_key(preferred) {
            return Some(preferred.to_owned());
        }
    }

    let mut best_channel = None;
    let mut best_score = f64::NEG_INFINITY;
    for channel in record.data_channels() {
        let Some(trace) = record.channel_values(&channel) else {
            continue;
        };
        let peaks = find_peaks(&trace, 100.0, 8);
        if peaks.len() < 12 {
            continue;
        }
        let score = quadratic_fit_r2(
            &(0..peaks.len()).map(|idx| idx as f64).collect::<Vec<_>>(),
            &peaks
                .iter()
                .map(|peak| peak.index as f64)
                .collect::<Vec<_>>(),
        );
        if score > best_score {
            best_score = score;
            best_channel = Some(channel);
        }
    }
    best_channel
}

fn quadratic_fit_r2(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.len() < 3 {
        return f64::NEG_INFINITY;
    }
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut s4 = 0.0;
    let mut t0 = 0.0;
    let mut t1 = 0.0;
    let mut t2 = 0.0;
    for (&xv, &yv) in x.iter().zip(y.iter()) {
        let x2 = xv * xv;
        s0 += 1.0;
        s1 += xv;
        s2 += x2;
        s3 += x2 * xv;
        s4 += x2 * x2;
        t0 += yv;
        t1 += xv * yv;
        t2 += x2 * yv;
    }

    let det = determinant3([[s4, s3, s2], [s3, s2, s1], [s2, s1, s0]]);
    if det.abs() <= f64::EPSILON {
        return f64::NEG_INFINITY;
    }
    let a = determinant3([[t2, s3, s2], [t1, s2, s1], [t0, s1, s0]]) / det;
    let b = determinant3([[s4, t2, s2], [s3, t1, s1], [s2, t0, s0]]) / det;
    let c = determinant3([[s4, s3, t2], [s3, s2, t1], [s2, s1, t0]]) / det;

    let mean_y = y.iter().sum::<f64>() / y.len() as f64;
    let mut ss_tot = 0.0;
    let mut ss_res = 0.0;
    for (&xv, &yv) in x.iter().zip(y.iter()) {
        let predicted = a * xv * xv + b * xv + c;
        let diff_tot = yv - mean_y;
        let diff_res = yv - predicted;
        ss_tot += diff_tot * diff_tot;
        ss_res += diff_res * diff_res;
    }
    if ss_tot <= f64::EPSILON {
        return f64::NEG_INFINITY;
    }
    1.0 - (ss_res / ss_tot)
}

fn determinant3(matrix: [[f64; 3]; 3]) -> f64 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

fn suggested_ladder_kind(
    record: &AbifRecord,
    size_standard_channel: &str,
    analysis_kind: Option<&AnalysisKind>,
) -> LadderKind {
    match analysis_kind {
        Some(AnalysisKind::Flt3) => LadderKind::Gs500Rox,
        Some(AnalysisKind::Clonality) if size_standard_channel == "DATA105" => {
            LadderKind::Liz500250
        }
        Some(AnalysisKind::Clonality) if size_standard_channel == "DATA4" => LadderKind::Rox400Hd,
        _ if record.tags.contains_key("DATA105") => LadderKind::Liz500250,
        _ => LadderKind::Rox400Hd,
    }
}

fn preferred_sample_channel(
    data_channels: &[String],
    size_standard_channel: &str,
    file_name: &str,
    analysis_kind: Option<&AnalysisKind>,
) -> Option<String> {
    if matches!(analysis_kind, Some(AnalysisKind::Flt3)) {
        if let Some(assay) = detect_flt3_assay(file_name) {
            return assay
                .channels
                .iter()
                .find_map(|channel| {
                    data_channels
                        .iter()
                        .find(|candidate| candidate.as_str() == *channel)
                        .cloned()
                })
                .or_else(|| {
                    data_channels
                        .iter()
                        .find(|channel| channel.as_str() != size_standard_channel)
                        .cloned()
                });
        }
    }
    None
}

fn select_ladder_peaks(
    raw_trace: &[f64],
    corrected_trace: &[f64],
    quantile_corrected_trace: &[f64],
    raw_min_height: f64,
    corrected_min_height: f64,
    min_distance: usize,
    max_peaks: usize,
    expected_peak_count: usize,
) -> Vec<Peak> {
    let target_candidate_count = expected_peak_count.min(8).max(4);
    let min_candidate_span = if expected_peak_count >= 20 { 1100 } else { 850 };
    let raw_candidates = adaptive_top_peak_candidates(
        raw_trace,
        raw_min_height,
        min_distance,
        max_peaks,
        target_candidate_count,
        min_candidate_span,
    );
    let corrected_candidates = adaptive_top_peak_candidates(
        corrected_trace,
        corrected_min_height,
        min_distance,
        max_peaks,
        target_candidate_count,
        min_candidate_span,
    );
    let quantile_candidates = adaptive_top_peak_candidates(
        quantile_corrected_trace,
        corrected_min_height,
        min_distance,
        max_peaks,
        target_candidate_count,
        min_candidate_span,
    );

    let merged_candidates = merge_candidate_sets(
        &[raw_candidates.clone(), corrected_candidates.clone(), quantile_candidates.clone()],
        min_distance,
        max_peaks,
    );
    let coverage_candidates = coverage_peak_candidates(
        quantile_corrected_trace,
        min_distance,
        max_peaks,
        expected_peak_count,
    );
    let merged_with_coverage = merge_candidate_sets(
        &[merged_candidates.clone(), coverage_candidates],
        min_distance,
        max_peaks,
    );

    let min_viable_raw = expected_peak_count.min(6).max(3);
    let raw_tail_ok = candidate_has_tail_coverage(&raw_candidates, expected_peak_count);
    if candidate_pool_rank(&merged_with_coverage) >= candidate_pool_rank(&merged_candidates) {
        return merged_with_coverage;
    }
    if merged_candidates.len() >= raw_candidates.len().max(min_viable_raw) && !raw_tail_ok {
        return merged_candidates;
    }
    if merged_candidates.len() >= raw_candidates.len().max(min_viable_raw) {
        return merged_candidates;
    }
    if raw_candidates.len() >= min_viable_raw && raw_tail_ok {
        return raw_candidates;
    }
    if corrected_candidates.len() >= quantile_candidates.len()
        && corrected_candidates.len() > raw_candidates.len()
    {
        return corrected_candidates;
    }
    if quantile_candidates.len() > raw_candidates.len() {
        return quantile_candidates;
    }
    raw_candidates
}

fn coverage_peak_candidates(
    values: &[f64],
    min_distance: usize,
    max_peaks: usize,
    expected_peak_count: usize,
) -> Vec<Peak> {
    if values.len() < 500 {
        return Vec::new();
    }

    let (window_start, window_end) = if expected_peak_count >= 20 {
        (1300usize, 4300usize)
    } else {
        (1150usize, 4300usize)
    };
    let hard_start = window_start.min(values.len().saturating_sub(1));
    let hard_end = window_end.min(values.len());
    if hard_end <= hard_start + 2 {
        return Vec::new();
    }

    let low_threshold = if expected_peak_count >= 20 { 12.0 } else { 10.0 };
    let all_peaks = find_peaks(values, low_threshold, min_distance);
    if all_peaks.is_empty() {
        return Vec::new();
    }

    let bucket_count = expected_peak_count.min(12).max(6);
    let bucket_width = ((hard_end - hard_start) as f64 / bucket_count as f64).ceil().max(1.0) as usize;
    let mut selected: Vec<Peak> = Vec::new();

    for bucket in 0..bucket_count {
        let start = hard_start.saturating_add(bucket * bucket_width);
        let end = if bucket + 1 >= bucket_count {
            hard_end
        } else {
            hard_start.saturating_add((bucket + 1) * bucket_width).min(hard_end)
        };
        let best = all_peaks
            .iter()
            .filter(|peak| peak.index >= start && peak.index < end)
            .max_by(|left, right| {
                left
                    .score
                    .partial_cmp(&right.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .cloned();
        if let Some(peak) = best {
            selected.push(peak);
        }
    }

    if selected.is_empty() {
        return Vec::new();
    }
    selected = select_diverse_peak_subset(selected, max_peaks.min(bucket_count + 4));
    selected.sort_by_key(|peak| peak.index);
    selected
}

fn adaptive_top_peak_candidates(
    values: &[f64],
    min_height: f64,
    min_distance: usize,
    max_peaks: usize,
    target_candidate_count: usize,
    min_candidate_span: usize,
) -> Vec<Peak> {
    let mut best = Vec::new();
    let mut thresholds = vec![
        min_height,
        (min_height * 0.80).max(20.0),
        (min_height * 0.60).max(20.0),
        (min_height * 0.45).max(20.0),
        20.0,
    ];
    thresholds.dedup_by(|left, right| (*left - *right).abs() < f64::EPSILON);

    for threshold in thresholds {
        let candidates = top_peak_candidates(values, threshold, min_distance, max_peaks);
        if candidate_pool_rank(&candidates) > candidate_pool_rank(&best) {
            best = candidates.clone();
        }
        if candidates.len() >= target_candidate_count && candidate_span(&candidates) >= min_candidate_span {
            return candidates;
        }
    }
    best
}

fn candidate_span(peaks: &[Peak]) -> usize {
    match (peaks.first(), peaks.last()) {
        (Some(first), Some(last)) => last.index.saturating_sub(first.index),
        _ => 0,
    }
}

fn candidate_pool_rank(peaks: &[Peak]) -> usize {
    peaks.len().saturating_mul(10_000).saturating_add(candidate_span(peaks))
}

fn candidate_has_tail_coverage(peaks: &[Peak], expected_peak_count: usize) -> bool {
    let Some(last) = peaks.last() else {
        return false;
    };
    if expected_peak_count >= 20 {
        return last.index >= 3300 && candidate_span(peaks) >= 1850;
    }
    last.index >= 4300 && candidate_span(peaks) >= 2400
}

fn top_peak_candidates(
    values: &[f64],
    min_height: f64,
    min_distance: usize,
    max_peaks: usize,
) -> Vec<Peak> {
    let mut peaks = find_peaks(values, min_height, min_distance);
    peaks = select_diverse_peak_subset(peaks, max_peaks);
    peaks.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .prominence
                    .partial_cmp(&left.prominence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                right
                    .height
                    .partial_cmp(&left.height)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.index.cmp(&right.index))
    });
    if peaks.len() > max_peaks {
        peaks.truncate(max_peaks);
    }
    peaks.sort_by_key(|peak| peak.index);
    peaks
}

fn select_diverse_peak_subset(peaks: Vec<Peak>, max_peaks: usize) -> Vec<Peak> {
    select_diverse_peak_subset_with_buckets(peaks, max_peaks, max_peaks.clamp(4, 8))
}

fn select_diverse_peak_subset_with_buckets(
    peaks: Vec<Peak>,
    max_peaks: usize,
    bucket_count: usize,
) -> Vec<Peak> {
    if peaks.len() <= max_peaks {
        return peaks;
    }

    let mut ranked = peaks;
    ranked.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .prominence
                    .partial_cmp(&left.prominence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                right
                    .height
                    .partial_cmp(&left.height)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.index.cmp(&right.index))
    });

    let min_index = ranked.iter().map(|peak| peak.index).min().unwrap_or(0);
    let max_index = ranked.iter().map(|peak| peak.index).max().unwrap_or(min_index);
    let span = max_index.saturating_sub(min_index);
    if span == 0 {
        ranked.truncate(max_peaks);
        return ranked;
    }

    let bucket_count = bucket_count.clamp(4, max_peaks.max(4));
    let bucket_width = ((span as f64) / bucket_count as f64).ceil().max(1.0) as usize;
    let mut selected: Vec<Peak> = Vec::new();

    for bucket in 0..bucket_count {
        let start = min_index.saturating_add(bucket * bucket_width);
        let end = if bucket + 1 >= bucket_count {
            max_index.saturating_add(1)
        } else {
            min_index.saturating_add((bucket + 1) * bucket_width)
        };
        if let Some(best) = ranked
            .iter()
            .filter(|peak| peak.index >= start && peak.index < end)
            .min_by(|left, right| {
                right
                    .score
                    .partial_cmp(&left.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .cloned()
        {
            if !selected.iter().any(|peak| peak.index == best.index) {
                selected.push(best);
            }
        }
    }

    for peak in ranked {
        if selected.len() >= max_peaks {
            break;
        }
        if selected.iter().any(|kept| kept.index == peak.index) {
            continue;
        }
        selected.push(peak);
    }

    selected
}

fn merge_candidate_sets(
    candidate_sets: &[Vec<Peak>],
    min_distance: usize,
    max_peaks: usize,
) -> Vec<Peak> {
    let mut combined = candidate_sets
        .iter()
        .flat_map(|peaks| peaks.iter().cloned())
        .collect::<Vec<_>>();
    combined.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .prominence
                    .partial_cmp(&left.prominence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                right
                    .height
                    .partial_cmp(&left.height)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.index.cmp(&right.index))
    });

    let merge_distance = ((min_distance as f64) * 0.25).floor() as usize;
    let mut merged: Vec<Peak> = Vec::new();
    'candidate: for candidate in combined {
        for kept in &merged {
            if candidate.index.abs_diff(kept.index) <= merge_distance {
                continue 'candidate;
            }
        }
        merged.push(candidate);
    }

    if candidate_span(&merged) < 1500 {
        merged = select_diverse_peak_subset(merged, max_peaks);
    } else if merged.len() > max_peaks {
        merged = select_diverse_peak_subset(merged, max_peaks);
    }
    merged.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .prominence
                    .partial_cmp(&left.prominence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                right
                    .height
                    .partial_cmp(&left.height)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.index.cmp(&right.index))
    });
    if merged.len() > max_peaks {
        merged.truncate(max_peaks);
    }
    merged.sort_by_key(|peak| peak.index);
    merged
}

fn build_ladder_fit_preview(
    ladder_peaks: &[Peak],
    sample_trace: &[f64],
    ladder: LadderKind,
) -> Option<LadderFitPreview> {
    let target_len = ladder.expected_peak_count();
    if ladder_peaks.len() < 2 || target_len < 2 {
        return None;
    }

    let peak_indices = ladder_peaks
        .iter()
        .map(|peak| peak.index)
        .collect::<Vec<_>>();
    let mut max_allowed_peak_gap = estimate_max_allowed_peak_gap(&peak_indices, 5.0);
    let mut gap_expansions = 0usize;
    let mut estimated_combination_count = 0usize;
    let mut candidate_generation_capped = false;
    let mut combinations = Vec::new();

    for expansion in 0..=LADDER_MAX_GAP_EXPANSIONS {
        estimated_combination_count = estimate_combination_count_capped(
            &peak_indices,
            target_len,
            max_allowed_peak_gap,
            MAX_CANDIDATE_COMBINATIONS + 1,
        );
        candidate_generation_capped = estimated_combination_count > MAX_CANDIDATE_COMBINATIONS;
        gap_expansions = expansion;
        if candidate_generation_capped {
            break;
        }

        combinations = if estimated_combination_count > BEAM_SEARCH_TRIGGER_COMBINATIONS {
            let peak_feature_by_index = ladder_peaks
                .iter()
                .map(|peak| (peak.index, peak.clone()))
                .collect::<BTreeMap<_, _>>();
            generate_peak_combinations_beam(
                &peak_indices,
                target_len,
                max_allowed_peak_gap,
                ladder.sizes(),
                &peak_feature_by_index,
                BEAM_SEARCH_WIDTH,
                BEAM_SEARCH_FINAL_CAP,
            )
        } else {
            generate_peak_combinations(
                &peak_indices,
                target_len,
                max_allowed_peak_gap,
                MAX_CANDIDATE_COMBINATIONS,
            )
        };
        if !combinations.is_empty() {
            break;
        }
        max_allowed_peak_gap = max_allowed_peak_gap.saturating_add(LADDER_GAP_EXPANSION_STEP);
    }

    if candidate_generation_capped {
        for keep_limit in [target_len + 8, target_len + 6, target_len + 4, target_len + 2] {
            for collapse_distance in [20usize, 28, 36, 48, 64, 80] {
                let reduced_ladder_peaks =
                    thin_peak_pool_for_ladder(ladder_peaks, target_len, collapse_distance, keep_limit);
                if reduced_ladder_peaks.len() < target_len
                    || reduced_ladder_peaks.len() >= ladder_peaks.len()
                {
                    continue;
                }
                let reduced_indices = reduced_ladder_peaks
                    .iter()
                    .map(|peak| peak.index)
                    .collect::<Vec<_>>();
                let reduced_max_gap = estimate_max_allowed_peak_gap(&reduced_indices, 4.0);
                let reduced_estimate = estimate_combination_count_capped(
                    &reduced_indices,
                    target_len,
                    reduced_max_gap,
                    MAX_CANDIDATE_COMBINATIONS + 1,
                );
                let peak_feature_by_index = reduced_ladder_peaks
                    .iter()
                    .map(|peak| (peak.index, peak.clone()))
                    .collect::<BTreeMap<_, _>>();
                let reduced_combinations =
                    if reduced_estimate > BEAM_SEARCH_TRIGGER_COMBINATIONS {
                        generate_peak_combinations_beam(
                            &reduced_indices,
                            target_len,
                            reduced_max_gap,
                            ladder.sizes(),
                            &peak_feature_by_index,
                            BEAM_SEARCH_WIDTH,
                            BEAM_SEARCH_FINAL_CAP,
                        )
                    } else if reduced_estimate <= MAX_CANDIDATE_COMBINATIONS {
                        generate_peak_combinations(
                            &reduced_indices,
                            target_len,
                            reduced_max_gap,
                            MAX_CANDIDATE_COMBINATIONS,
                        )
                    } else {
                        continue;
                    };
                if reduced_combinations.is_empty() {
                    continue;
                }
                let ladder_sizes = ladder.sizes();
                let mut best = select_best_combination(
                    &reduced_combinations,
                    ladder_sizes,
                    ladder,
                    &reduced_ladder_peaks,
                );
                let mut sizing_model = best.as_ref().and_then(|entry| {
                    fit_best_sizing_model(&entry.indices, ladder_sizes, sample_trace)
                });
                let mut refinement = None;

                if let (Some(best_entry), Some(model)) = (best.as_ref(), sizing_model.as_ref()) {
                    if let Some(refined) = refine_best_combination(
                        &reduced_indices,
                        &best_entry.indices,
                        ladder_sizes,
                        model,
                    ) {
                        let refined_qr2 = quadratic_fit_r2(
                            ladder_sizes,
                            &refined
                                .refined_scan_indices
                                .iter()
                                .map(|value| *value as f64)
                                .collect::<Vec<_>>(),
                        );
                        let refined_curvature =
                            curvature_score(ladder_sizes, &refined.refined_scan_indices);
                        let refined_domain_penalty = ladder_domain_penalty(
                            ladder,
                            &refined.refined_scan_indices,
                            &peak_feature_by_index,
                            &reduced_ladder_peaks,
                        );
                        let refined_peak_penalty = ladder_peak_sequence_penalty(
                            &refined.refined_scan_indices,
                            ladder_sizes,
                            &peak_feature_by_index,
                            &reduced_ladder_peaks,
                        );
                        best = Some(CombinationScore {
                            indices: refined.refined_scan_indices.clone(),
                            curvature_score: refined_curvature,
                            quadratic_r2: refined_qr2,
                            domain_penalty: refined_domain_penalty,
                            peak_penalty: refined_peak_penalty,
                            blended_score: refined_curvature
                                + refined_domain_penalty
                                + refined_peak_penalty,
                        });
                        sizing_model = fit_best_sizing_model(
                            &refined.refined_scan_indices,
                            ladder_sizes,
                            sample_trace,
                        );
                        refinement = Some(RefinementPreview {
                            changed_step_indices: refined.changed_step_indices,
                            original_scan_indices: refined.original_scan_indices,
                            refined_scan_indices: refined.refined_scan_indices,
                            refined_curvature_score: refined_curvature,
                            refined_quadratic_r2: refined_qr2,
                        });
                    }
                }

                return Some(LadderFitPreview {
                    max_allowed_peak_gap: reduced_max_gap,
                    gap_expansions,
                    estimated_combination_count: reduced_estimate,
                    candidate_generation_capped: false,
                    evaluated_combination_count: reduced_combinations.len(),
                    best_scan_indices: best
                        .as_ref()
                        .map(|entry| entry.indices.clone())
                        .unwrap_or_default(),
                    best_curvature_score: best.as_ref().map(|entry| entry.curvature_score),
                    best_quadratic_r2: best.as_ref().map(|entry| entry.quadratic_r2),
                    sizing_model,
                    refinement,
                });
            }
        }
        return Some(LadderFitPreview {
            max_allowed_peak_gap,
            gap_expansions,
            estimated_combination_count,
            candidate_generation_capped: true,
            evaluated_combination_count: 0,
            best_scan_indices: Vec::new(),
            best_curvature_score: None,
            best_quadratic_r2: None,
            sizing_model: None,
            refinement: None,
        });
    }
    let ladder_sizes = ladder.sizes();
    let peak_feature_by_index = ladder_peaks
        .iter()
        .map(|peak| (peak.index, peak.clone()))
        .collect::<BTreeMap<_, _>>();
    let mut best = select_best_combination(&combinations, ladder_sizes, ladder, ladder_peaks);
    let mut sizing_model = best
        .as_ref()
        .and_then(|entry| fit_best_sizing_model(&entry.indices, ladder_sizes, sample_trace));
    let mut refinement = None;

    if let (Some(best_entry), Some(model)) = (best.as_ref(), sizing_model.as_ref()) {
        if let Some(refined) =
            refine_best_combination(&peak_indices, &best_entry.indices, ladder_sizes, model)
        {
            let refined_qr2 = quadratic_fit_r2(
                ladder_sizes,
                &refined
                    .refined_scan_indices
                    .iter()
                    .map(|value| *value as f64)
                    .collect::<Vec<_>>(),
            );
            let refined_curvature = curvature_score(ladder_sizes, &refined.refined_scan_indices);
            let refined_domain_penalty =
                ladder_domain_penalty(
                    ladder,
                    &refined.refined_scan_indices,
                    &peak_feature_by_index,
                    ladder_peaks,
                );
            let refined_peak_penalty = ladder_peak_sequence_penalty(
                &refined.refined_scan_indices,
                ladder_sizes,
                &peak_feature_by_index,
                ladder_peaks,
            );
            best = Some(CombinationScore {
                indices: refined.refined_scan_indices.clone(),
                curvature_score: refined_curvature,
                quadratic_r2: refined_qr2,
                domain_penalty: refined_domain_penalty,
                peak_penalty: refined_peak_penalty,
                blended_score: refined_curvature + refined_domain_penalty + refined_peak_penalty,
            });
            sizing_model =
                fit_best_sizing_model(&refined.refined_scan_indices, ladder_sizes, sample_trace);
            refinement = Some(RefinementPreview {
                changed_step_indices: refined.changed_step_indices,
                original_scan_indices: refined.original_scan_indices,
                refined_scan_indices: refined.refined_scan_indices,
                refined_curvature_score: refined_curvature,
                refined_quadratic_r2: refined_qr2,
            });
        }
    }

    Some(LadderFitPreview {
        max_allowed_peak_gap,
        gap_expansions,
        estimated_combination_count,
        candidate_generation_capped: false,
        evaluated_combination_count: combinations.len(),
        best_scan_indices: best
            .as_ref()
            .map(|entry| entry.indices.clone())
            .unwrap_or_default(),
        best_curvature_score: best.as_ref().map(|entry| entry.curvature_score),
        best_quadratic_r2: best.as_ref().map(|entry| entry.quadratic_r2),
        sizing_model,
        refinement,
    })
}

fn thin_peak_pool_for_ladder(
    ladder_peaks: &[Peak],
    target_len: usize,
    collapse_distance: usize,
    keep_limit: usize,
) -> Vec<Peak> {
    if ladder_peaks.len() <= keep_limit {
        return ladder_peaks.to_vec();
    }

    let diverse = select_diverse_peak_subset_with_buckets(
        ladder_peaks.to_vec(),
        keep_limit,
        keep_limit.clamp(8, 24),
    );
    if diverse.len() < target_len {
        return ladder_peaks.to_vec();
    }

    let mut ranked = diverse;
    ranked.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .prominence
                    .partial_cmp(&left.prominence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                right
                    .height
                    .partial_cmp(&left.height)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.index.cmp(&right.index))
    });

    let mut reduced: Vec<Peak> = Vec::new();
    'candidate: for peak in ranked {
        for kept in &reduced {
            if peak.index.abs_diff(kept.index) <= collapse_distance {
                continue 'candidate;
            }
        }
        reduced.push(peak);
        if reduced.len() >= keep_limit {
            break;
        }
    }

    if reduced.len() < target_len {
        return ladder_peaks.to_vec();
    }
    reduced.sort_by_key(|peak| peak.index);
    reduced
}

fn estimate_max_allowed_peak_gap(peak_indices: &[usize], multiplier: f64) -> usize {
    if peak_indices.len() < 2 {
        return 0;
    }
    let mean_gap = peak_indices
        .windows(2)
        .map(|window| window[1].saturating_sub(window[0]) as f64)
        .sum::<f64>()
        / (peak_indices.len() - 1) as f64;
    (mean_gap * multiplier).round().max(1.0) as usize
}

fn estimate_combination_count_capped(
    peak_indices: &[usize],
    target_len: usize,
    max_gap: usize,
    cap: usize,
) -> usize {
    fn dfs(
        peak_indices: &[usize],
        start: usize,
        chosen: usize,
        target_len: usize,
        max_gap: usize,
        last_selected: Option<usize>,
        cap: usize,
        memo: &mut BTreeMap<(usize, usize, Option<usize>), usize>,
    ) -> usize {
        if chosen == target_len {
            return 1;
        }
        if start >= peak_indices.len() {
            return 0;
        }
        let key = (start, chosen, last_selected);
        if let Some(value) = memo.get(&key) {
            return *value;
        }

        let mut total = 0usize;
        for candidate_index in start..peak_indices.len() {
            if let Some(last_index) = last_selected {
                let gap = peak_indices[candidate_index].abs_diff(peak_indices[last_index]);
                if gap > max_gap {
                    continue;
                }
            }
            total = total.saturating_add(dfs(
                peak_indices,
                candidate_index + 1,
                chosen + 1,
                target_len,
                max_gap,
                Some(candidate_index),
                cap,
                memo,
            ));
            if total >= cap {
                total = cap;
                break;
            }
        }

        memo.insert(key, total);
        total
    }

    if target_len == 0 || peak_indices.len() < target_len {
        return 0;
    }
    let mut memo = BTreeMap::new();
    dfs(
        peak_indices,
        0,
        0,
        target_len,
        max_gap,
        None,
        cap,
        &mut memo,
    )
}

fn generate_peak_combinations(
    peak_indices: &[usize],
    target_len: usize,
    max_gap: usize,
    max_combinations: usize,
) -> Vec<Vec<usize>> {
    fn dfs(
        peak_indices: &[usize],
        start: usize,
        target_len: usize,
        max_gap: usize,
        current: &mut Vec<usize>,
        results: &mut Vec<Vec<usize>>,
        max_combinations: usize,
    ) {
        if results.len() >= max_combinations {
            return;
        }
        if current.len() == target_len {
            results.push(current.clone());
            return;
        }
        if start >= peak_indices.len() {
            return;
        }

        for candidate_index in start..peak_indices.len() {
            if let Some(&last_scan) = current.last() {
                let gap = peak_indices[candidate_index].abs_diff(last_scan);
                if gap > max_gap {
                    continue;
                }
            }
            current.push(peak_indices[candidate_index]);
            dfs(
                peak_indices,
                candidate_index + 1,
                target_len,
                max_gap,
                current,
                results,
                max_combinations,
            );
            current.pop();
            if results.len() >= max_combinations {
                return;
            }
        }
    }

    if target_len == 0 || peak_indices.len() < target_len {
        return Vec::new();
    }

    let mut results = Vec::new();
    let mut current = Vec::with_capacity(target_len);
    dfs(
        peak_indices,
        0,
        target_len,
        max_gap,
        &mut current,
        &mut results,
        max_combinations,
    );
    results
}

fn generate_peak_combinations_beam(
    peak_indices: &[usize],
    target_len: usize,
    max_gap: usize,
    ladder_sizes: &[f64],
    peak_feature_by_index: &BTreeMap<usize, Peak>,
    beam_width: usize,
    final_cap: usize,
) -> Vec<Vec<usize>> {
    #[derive(Clone)]
    struct BeamState {
        indices: Vec<usize>,
        next_start: usize,
        score: f64,
    }

    if target_len == 0 || peak_indices.len() < target_len || ladder_sizes.len() < target_len {
        return Vec::new();
    }

    let mut states = vec![BeamState {
        indices: Vec::with_capacity(target_len),
        next_start: 0,
        score: 0.0,
    }];

    for step in 0..target_len {
        let remaining_after_pick = target_len.saturating_sub(step + 1);
        let mut next_states = Vec::new();

        for state in &states {
            if peak_indices.len().saturating_sub(state.next_start) < remaining_after_pick + 1 {
                continue;
            }

            for candidate_index in state.next_start..peak_indices.len() {
                if peak_indices.len().saturating_sub(candidate_index + 1) < remaining_after_pick {
                    break;
                }
                if let Some(&last_scan) = state.indices.last() {
                    let gap = peak_indices[candidate_index].abs_diff(last_scan);
                    if gap > max_gap {
                        continue;
                    }
                }

                let mut indices = state.indices.clone();
                indices.push(peak_indices[candidate_index]);
                let score = partial_combination_beam_score(
                    &indices,
                    &ladder_sizes[..indices.len()],
                    peak_feature_by_index,
                );
                next_states.push(BeamState {
                    indices,
                    next_start: candidate_index + 1,
                    score,
                });
            }
        }

        if next_states.is_empty() {
            return Vec::new();
        }

        next_states.sort_by(|left, right| {
            left.score
                .partial_cmp(&right.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| left.indices.cmp(&right.indices))
        });
        if next_states.len() > beam_width {
            next_states.truncate(beam_width);
        }
        states = next_states;
    }

    let mut results = states
        .into_iter()
        .map(|state| (state.score, state.indices))
        .collect::<Vec<_>>();
    results.sort_by(|left, right| {
        left.0
            .partial_cmp(&right.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.1.cmp(&right.1))
    });
    if results.len() > final_cap {
        results.truncate(final_cap);
    }
    results.into_iter().map(|(_, indices)| indices).collect()
}

fn partial_combination_beam_score(
    scans: &[usize],
    ladder_sizes: &[f64],
    peak_feature_by_index: &BTreeMap<usize, Peak>,
) -> f64 {
    if scans.is_empty() || ladder_sizes.is_empty() || scans.len() != ladder_sizes.len() {
        return f64::INFINITY;
    }

    let prefix_curvature = if scans.len() >= 3 {
        curvature_score(ladder_sizes, scans)
    } else {
        0.0
    };

    let gap_ratios = scans
        .windows(2)
        .zip(ladder_sizes.windows(2))
        .filter_map(|(scan_window, bp_window)| {
            let bp_gap = bp_window[1] - bp_window[0];
            if bp_gap <= f64::EPSILON {
                None
            } else {
                Some((scan_window[1] as f64 - scan_window[0] as f64) / bp_gap)
            }
        })
        .collect::<Vec<_>>();
    let gap_penalty = coefficient_of_variation_penalty(&gap_ratios, 0.45) * 0.25;

    let score_reward = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.score.max(1.0).ln()))
        .sum::<f64>()
        / scans.len() as f64;

    let purity_penalty = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan))
        .map(|peak| {
            let height = peak.height.max(1.0);
            let baseline_ratio = (peak.local_baseline.max(0.0) / height).clamp(0.0, 1.5);
            let purity = (peak.prominence / height).clamp(0.0, 1.0);
            baseline_ratio * (1.0 - 0.60 * purity)
        })
        .sum::<f64>()
        / scans.len() as f64;

    let early_dense_penalty = if gap_ratios.len() >= 2 {
        let early_len = gap_ratios.len().min(3);
        let median_gap_ratio = median(&gap_ratios).max(1.0);
        gap_ratios
            .iter()
            .take(early_len)
            .map(|ratio| ((median_gap_ratio * 0.68 - *ratio).max(0.0)) / median_gap_ratio)
            .sum::<f64>()
            / early_len as f64
    } else {
        0.0
    };

    prefix_curvature + gap_penalty + purity_penalty * 0.60 + early_dense_penalty * 0.70
        - score_reward * 0.04
}

#[derive(Debug, Clone, PartialEq)]
struct CombinationScore {
    indices: Vec<usize>,
    curvature_score: f64,
    quadratic_r2: f64,
    domain_penalty: f64,
    peak_penalty: f64,
    blended_score: f64,
}

fn select_best_combination(
    combinations: &[Vec<usize>],
    ladder_sizes: &[f64],
    ladder: LadderKind,
    peak_features: &[Peak],
) -> Option<CombinationScore> {
    let peak_feature_by_index = peak_features
        .iter()
        .map(|peak| (peak.index, peak.clone()))
        .collect::<BTreeMap<_, _>>();
    combinations
        .iter()
        .filter(|combo| combo.len() == ladder_sizes.len())
        .map(|combo| {
            let curvature = curvature_score(ladder_sizes, combo);
            let qr2 = quadratic_fit_r2(
                ladder_sizes,
                &combo.iter().map(|value| *value as f64).collect::<Vec<_>>(),
            );
            let domain_penalty =
                ladder_domain_penalty(ladder, combo, &peak_feature_by_index, peak_features);
            let peak_penalty =
                ladder_peak_sequence_penalty(combo, ladder_sizes, &peak_feature_by_index, peak_features);
            CombinationScore {
                indices: combo.clone(),
                curvature_score: curvature,
                quadratic_r2: qr2,
                domain_penalty,
                peak_penalty,
                blended_score: curvature + domain_penalty + peak_penalty,
            }
        })
        .min_by(|left, right| {
            left.blended_score
                .partial_cmp(&right.blended_score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    left.peak_penalty
                        .partial_cmp(&right.peak_penalty)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    left.curvature_score
                        .partial_cmp(&right.curvature_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    right
                        .quadratic_r2
                        .partial_cmp(&left.quadratic_r2)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| left.indices.cmp(&right.indices))
        })
}

fn ladder_domain_penalty(
    ladder: LadderKind,
    scans: &[usize],
    peak_feature_by_index: &BTreeMap<usize, Peak>,
    all_peak_features: &[Peak],
) -> f64 {
    if scans.is_empty() {
        return f64::INFINITY;
    }

    let (hard_min, hard_max, preferred_min, preferred_max, max_first_anchor) = match ladder {
        LadderKind::Rox400Hd => (
            ROX_HARD_TIME_MIN,
            ROX_HARD_TIME_MAX,
            ROX_PREFERRED_TIME_MIN,
            ROX_PREFERRED_TIME_MAX,
            ROX_MAX_FIRST_ANCHOR,
        ),
        LadderKind::Gs500Rox => (
            GS500ROX_HARD_TIME_MIN,
            GS500ROX_HARD_TIME_MAX,
            GS500ROX_PREFERRED_TIME_MIN,
            GS500ROX_PREFERRED_TIME_MAX,
            GS500ROX_MAX_FIRST_ANCHOR,
        ),
        LadderKind::Liz500250 => (
            LIZ_HARD_TIME_MIN,
            LIZ_HARD_TIME_MAX,
            LIZ_PREFERRED_TIME_MIN,
            LIZ_PREFERRED_TIME_MAX,
            LIZ_MAX_FIRST_ANCHOR,
        ),
    };

    let scans_f = scans.iter().map(|value| *value as f64).collect::<Vec<_>>();
    let hard_out_fraction = scans_f
        .iter()
        .filter(|value| **value < hard_min || **value > hard_max)
        .count() as f64
        / scans_f.len() as f64;

    let median_scan = scans_f[scans_f.len() / 2];
    let median_window_penalty = if median_scan < preferred_min {
        (preferred_min - median_scan) / preferred_min.max(1.0)
    } else if median_scan > preferred_max {
        (median_scan - preferred_max) / preferred_max.max(1.0)
    } else {
        0.0
    };

    let first_anchor_late_penalty = if scans_f[0] > max_first_anchor {
        (scans_f[0] - max_first_anchor) / max_first_anchor.max(1.0)
    } else {
        0.0
    };
    let first_anchor_early_penalty = if scans_f[0] < preferred_min {
        (preferred_min - scans_f[0]) / preferred_min.max(1.0)
    } else {
        0.0
    };

    let gap_cv_penalty = if scans_f.len() >= 3 {
        let gaps = scans_f.windows(2).map(|window| window[1] - window[0]).collect::<Vec<_>>();
        let mean_gap = gaps.iter().sum::<f64>() / gaps.len() as f64;
        if mean_gap <= f64::EPSILON {
            1.0
        } else {
            let variance = gaps
                .iter()
                .map(|gap| {
                    let d = *gap - mean_gap;
                    d * d
                })
                .sum::<f64>()
                / gaps.len() as f64;
            let cv = variance.sqrt() / mean_gap;
            (cv - 0.55).max(0.0)
        }
    } else {
        0.0
    };

    let intensities = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.height))
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect::<Vec<_>>();
    let intensity_cv_penalty = if intensities.len() >= 4 {
        let mean_intensity = intensities.iter().sum::<f64>() / intensities.len() as f64;
        if mean_intensity <= f64::EPSILON {
            1.0
        } else {
            let variance = intensities
                .iter()
                .map(|value| {
                    let d = *value - mean_intensity;
                    d * d
                })
                .sum::<f64>()
                / intensities.len() as f64;
            let cv = variance.sqrt() / mean_intensity;
            (cv - 0.60).max(0.0)
        }
    } else {
        0.0
    };

    let prominences = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.prominence))
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect::<Vec<_>>();
    let widths = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.width))
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect::<Vec<_>>();

    let weak_prominence_penalty = if !prominences.is_empty() {
        let median_prominence = median(&prominences).max(1.0);
        let target_prominence = (median_prominence * 0.45).max(25.0);
        prominences
            .iter()
            .map(|value| ((target_prominence - *value).max(0.0)) / target_prominence)
            .sum::<f64>()
            / prominences.len() as f64
    } else {
        0.0
    };

    let width_cv_penalty = coefficient_of_variation_penalty(&widths, 0.85);

    let early_skip_penalty = if !all_peak_features.is_empty() {
        let early_peaks = all_peak_features
            .iter()
            .filter(|peak| {
                let scan = peak.index as f64;
                scan >= preferred_min - 50.0
                    && scan <= (max_first_anchor + 140.0)
                    && peak.prominence >= 30.0
            })
            .collect::<Vec<_>>();
        if early_peaks.is_empty() {
            0.0
        } else {
            let early_scores = early_peaks.iter().map(|peak| peak.score).collect::<Vec<_>>();
            let score_floor = median(&early_scores) * 0.65;
            let earliest_strong = early_peaks
                .iter()
                .filter(|peak| peak.score >= score_floor.max(20.0))
                .map(|peak| peak.index as f64)
                .min_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
            if let Some(earliest) = earliest_strong {
                ((scans_f[0] - earliest - 25.0).max(0.0)) / 150.0
            } else {
                0.0
            }
        }
    } else {
        0.0
    };

    let early_cluster_penalty = if scans_f.len() >= 3 {
        let early_scan_count = scans_f
            .iter()
            .take_while(|scan| **scan < preferred_min + 250.0)
            .count();
        if early_scan_count >= 2 {
            let early_gaps = scans_f[..early_scan_count]
                .windows(2)
                .map(|window| window[1] - window[0])
                .collect::<Vec<_>>();
            let all_gaps = scans_f
                .windows(2)
                .map(|window| window[1] - window[0])
                .collect::<Vec<_>>();
            if early_gaps.is_empty() || all_gaps.is_empty() {
                0.0
            } else {
                let early_gap_median = median(&early_gaps);
                let global_gap_median = median(&all_gaps).max(1.0);
                ((global_gap_median * 0.60 - early_gap_median).max(0.0)) / global_gap_median
            }
        } else {
            0.0
        }
    } else {
        0.0
    };

    let pre_window_penalty = scans_f
        .iter()
        .filter(|scan| **scan < preferred_min)
        .count() as f64
        / scans_f.len().max(1) as f64;

    hard_out_fraction * LADDER_DOMAIN_TIME_WEIGHT
        + median_window_penalty * LADDER_DOMAIN_TIME_WEIGHT
        + first_anchor_late_penalty * LADDER_DOMAIN_FIRST_ANCHOR_WEIGHT
        + first_anchor_early_penalty * (LADDER_DOMAIN_FIRST_ANCHOR_WEIGHT + 0.75)
        + gap_cv_penalty * LADDER_DOMAIN_GAP_WEIGHT
        + intensity_cv_penalty * LADDER_DOMAIN_INTENSITY_WEIGHT
        + weak_prominence_penalty * 0.35
        + width_cv_penalty * 0.20
        + early_skip_penalty * LADDER_DOMAIN_FIRST_ANCHOR_WEIGHT
        + early_cluster_penalty * 0.90
        + pre_window_penalty * 4.50
}

fn ladder_peak_sequence_penalty(
    scans: &[usize],
    ladder_sizes: &[f64],
    peak_feature_by_index: &BTreeMap<usize, Peak>,
    all_peak_features: &[Peak],
) -> f64 {
    if scans.len() < 2 || scans.len() != ladder_sizes.len() {
        return 0.0;
    }

    let gap_ratios = scans
        .windows(2)
        .zip(ladder_sizes.windows(2))
        .filter_map(|(scan_window, bp_window)| {
            let bp_gap = bp_window[1] - bp_window[0];
            if bp_gap <= f64::EPSILON {
                return None;
            }
            Some((scan_window[1] as f64 - scan_window[0] as f64) / bp_gap)
        })
        .collect::<Vec<_>>();
    let gap_ratio_penalty = coefficient_of_variation_penalty(&gap_ratios, 0.35) * 0.35;
    let median_gap_ratio = median(&gap_ratios).max(1.0);

    let gs500_start_penalty = if ladder_sizes.len() == 16 {
        let first_scan = scans[0] as f64;
        let first_anchor_penalty =
            ((GS500ROX_PREFERRED_TIME_MIN + 30.0 - first_scan).max(0.0)) / 95.0;
        let first_gap_penalty = scans
            .windows(2)
            .next()
            .map(|window| {
                let gap = window[1] as f64 - window[0] as f64;
                let low_penalty = ((66.0 - gap).max(0.0)) / 28.0;
                let high_penalty = ((gap - 80.0).max(0.0)) / 26.0;
                low_penalty + high_penalty
            })
            .unwrap_or(0.0);
        first_anchor_penalty * 1.35 + first_gap_penalty * 0.45
    } else {
        0.0
    };

    let bridge_skip_penalty = if ladder_sizes.len() == 16 && scans.len() >= 2 {
        let first = scans[0];
        let second = scans[1];
        let candidate = all_peak_features
            .iter()
            .filter(|peak| {
                peak.index > first.saturating_add(10)
                    && peak.index + 10 < second
                    && peak.prominence >= 120.0
                    && peak.score >= 850.0
            })
            .max_by(|left, right| {
                left.score
                    .partial_cmp(&right.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        if let Some(peak) = candidate {
            let score_ref = peak_feature_by_index
                .get(&scans[1])
                .map(|selected| selected.score.max(1.0))
                .unwrap_or(1.0);
            ((peak.score / score_ref).clamp(0.0, 1.25)) * 0.85
        } else {
            0.0
        }
    } else {
        0.0
    };

    let gs500_blob_cluster_penalty = if ladder_sizes.len() == 16 && scans.len() >= 4 {
        let cluster_candidates = all_peak_features
            .iter()
            .filter(|peak| {
                peak.index >= scans[0].saturating_sub(20)
                    && peak.index <= scans[2].saturating_add(20)
                    && peak.prominence >= 500.0
            })
            .collect::<Vec<_>>();
        let dense_cluster = cluster_candidates.len() >= 4;
        let early_start = scans[0] < 1650;
        if dense_cluster && early_start {
            ((1650.0 - scans[0] as f64).max(0.0) / 140.0) * 1.10
        } else {
            0.0
        }
    } else {
        0.0
    };

    let gs500_tail_twin_penalty = if ladder_sizes.len() == 16 && scans.len() >= 3 {
        let selected_tail = &scans[scans.len() - 3..];
        let strong_tail_candidates = all_peak_features
            .iter()
            .filter(|peak| peak.index >= selected_tail[0].saturating_sub(10) && peak.prominence >= 300.0)
            .map(|peak| peak.index)
            .collect::<Vec<_>>();
        let has_late_twin = strong_tail_candidates.iter().any(|index| *index >= scans[scans.len() - 1] + 35);
        let second_last_too_early = scans[scans.len() - 2] + 55 < scans[scans.len() - 1];
        if has_late_twin && second_last_too_early {
            0.95
        } else {
            0.0
        }
    } else {
        0.0
    };

    let early_gap_penalty = if gap_ratios.len() >= 3 {
        let early_len = gap_ratios.len().min(4);
        gap_ratios
            .iter()
            .take(early_len)
            .map(|ratio| ((median_gap_ratio * 0.72 - *ratio).max(0.0)) / median_gap_ratio)
            .sum::<f64>()
            / early_len as f64
            * 0.95
    } else {
        0.0
    };

    let selected_scores = scans
        .iter()
        .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.score))
        .collect::<Vec<_>>();
    let weak_score_penalty = if !selected_scores.is_empty() {
        let median_score = median(&selected_scores).max(1.0);
        let target_score = (median_score * 0.45).max(20.0);
        selected_scores
            .iter()
            .map(|value| ((target_score - *value).max(0.0)) / target_score)
            .sum::<f64>()
            / selected_scores.len() as f64
            * 0.55
    } else {
        0.0
    };

    let neighbor_intensity_penalty = scans
        .windows(2)
        .filter_map(|window| {
            let left = peak_feature_by_index.get(&window[0])?;
            let right = peak_feature_by_index.get(&window[1])?;
            let high = left.prominence.max(right.prominence).max(1.0);
            let low = left.prominence.min(right.prominence).max(1.0);
            Some((high / low).ln().max(0.0))
        })
        .sum::<f64>()
        / (scans.len() - 1) as f64
        * 0.20;

    let late_weak_penalty = if !selected_scores.is_empty() {
        let median_score = median(&selected_scores).max(1.0);
        let selected_heights = scans
            .iter()
            .filter_map(|scan| peak_feature_by_index.get(scan).map(|peak| peak.height))
            .collect::<Vec<_>>();
        let median_height = median(&selected_heights).max(1.0);
        let late_score_floor = (median_score * 0.30).max(16.0);
        let late_height_floor = (median_height * 0.35).max(55.0);
        scans.iter()
            .enumerate()
            .filter_map(|(index, scan)| {
                let peak = peak_feature_by_index.get(scan)?;
                let position = index as f64 / (scans.len().saturating_sub(1).max(1) as f64);
                if position < 0.55 {
                    return None;
                }
                let score_penalty = ((late_score_floor - peak.score).max(0.0)) / late_score_floor;
                let height_penalty =
                    ((late_height_floor - peak.height).max(0.0)) / late_height_floor;
                Some((score_penalty + height_penalty) * position)
            })
            .sum::<f64>()
            / scans.len() as f64
            * 0.70
    } else {
        0.0
    };

    let early_baseline_penalty = if scans.len() >= 4 {
        let early_len = scans.len().min(5);
        scans.iter()
            .take(early_len)
            .filter_map(|scan| {
                let peak = peak_feature_by_index.get(scan)?;
                let height = peak.height.max(1.0);
                let baseline_ratio = (peak.local_baseline.max(0.0) / height).clamp(0.0, 1.5);
                let purity = (peak.prominence / height).clamp(0.0, 1.0);
                Some((baseline_ratio * (1.0 - 0.65 * purity)).max(0.0))
            })
            .sum::<f64>()
            / early_len as f64
            * 0.85
    } else {
        0.0
    };

    gap_ratio_penalty
        + gs500_start_penalty
        + bridge_skip_penalty
        + gs500_blob_cluster_penalty
        + gs500_tail_twin_penalty
        + early_gap_penalty
        + weak_score_penalty
        + neighbor_intensity_penalty
        + late_weak_penalty
        + early_baseline_penalty
}

fn coefficient_of_variation_penalty(values: &[f64], tolerance: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    if mean <= f64::EPSILON {
        return 1.0;
    }
    let variance = values
        .iter()
        .map(|value| {
            let delta = *value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    let cv = variance.sqrt() / mean;
    (cv - tolerance).max(0.0)
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        0.5 * (sorted[mid - 1] + sorted[mid])
    } else {
        sorted[mid]
    }
}

fn curvature_score(ladder_sizes: &[f64], scans: &[usize]) -> f64 {
    if ladder_sizes.len() != scans.len() || ladder_sizes.len() < 3 {
        return f64::INFINITY;
    }

    ladder_sizes
        .windows(3)
        .zip(scans.windows(3))
        .map(|(x_window, y_window)| {
            let x0 = x_window[0];
            let x1 = x_window[1];
            let x2 = x_window[2];
            let y0 = y_window[0] as f64;
            let y1 = y_window[1] as f64;
            let y2 = y_window[2] as f64;
            let left_dx = x1 - x0;
            let right_dx = x2 - x1;
            let span = x2 - x0;
            if left_dx <= f64::EPSILON || right_dx <= f64::EPSILON || span <= f64::EPSILON {
                return f64::INFINITY;
            }
            let left_slope = (y1 - y0) / left_dx;
            let right_slope = (y2 - y1) / right_dx;
            (2.0 * (right_slope - left_slope) / span).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn fit_best_sizing_model(
    scans: &[usize],
    ladder_sizes: &[f64],
    sample_trace: &[f64],
) -> Option<SizingModelPreview> {
    if scans.len() != ladder_sizes.len() || scans.len() < 2 {
        return None;
    }

    let x = scans.iter().map(|value| *value as f64).collect::<Vec<_>>();
    if let Some(model) = fit_monotone_spline_sizing_model(&x, ladder_sizes, sample_trace) {
        return Some(model);
    }

    let max_degree = (scans.len().saturating_sub(1)).min(3);
    for degree in (1..=max_degree).rev() {
        let Some(coefficients) = fit_polynomial_least_squares(&x, ladder_sizes, degree) else {
            continue;
        };
        let predicted = x
            .iter()
            .map(|value| eval_polynomial(&coefficients, *value))
            .collect::<Vec<_>>();
        let qc_metrics = compute_ladder_qc_metrics(ladder_sizes, &predicted);
        if qc_metrics.monotonic_on_ladder {
            let sample_basepairs = (0..sample_trace.len())
                .map(|time| eval_polynomial(&coefficients, time as f64))
                .collect::<Vec<_>>();
            let sample_mapping = build_sample_mapping_preview(sample_trace, &sample_basepairs);
            return Some(SizingModelPreview {
                strategy: "polynomial_fallback".to_owned(),
                degree,
                coefficients,
                predicted_ladder_basepairs: predicted,
                qc_metrics,
                sample_mapping,
            });
        }
    }
    None
}

fn fit_monotone_spline_sizing_model(
    x: &[f64],
    ladder_sizes: &[f64],
    sample_trace: &[f64],
) -> Option<SizingModelPreview> {
    let tangents = monotone_cubic_tangents(x, ladder_sizes)?;
    let predicted = x
        .iter()
        .map(|value| eval_monotone_cubic_spline(x, ladder_sizes, &tangents, *value))
        .collect::<Vec<_>>();
    let qc_metrics = compute_ladder_qc_metrics(ladder_sizes, &predicted);
    if !qc_metrics.monotonic_on_ladder {
        return None;
    }

    let sample_basepairs = (0..sample_trace.len())
        .map(|time| eval_monotone_cubic_spline(x, ladder_sizes, &tangents, time as f64))
        .collect::<Vec<_>>();
    let sample_mapping = build_sample_mapping_preview(sample_trace, &sample_basepairs);

    Some(SizingModelPreview {
        strategy: "willros_monotone_spline".to_owned(),
        degree: 3,
        coefficients: tangents,
        predicted_ladder_basepairs: predicted,
        qc_metrics,
        sample_mapping,
    })
}

fn fit_polynomial_least_squares(x: &[f64], y: &[f64], degree: usize) -> Option<Vec<f64>> {
    if x.len() != y.len() || x.is_empty() {
        return None;
    }
    let order = degree + 1;
    let mut normal = vec![vec![0.0; order]; order];
    let mut rhs = vec![0.0; order];

    for row in 0..order {
        for col in 0..order {
            let power = (row + col) as i32;
            normal[row][col] = x.iter().map(|value| value.powi(power)).sum::<f64>();
        }
        rhs[row] = x
            .iter()
            .zip(y.iter())
            .map(|(x_value, y_value)| y_value * x_value.powi(row as i32))
            .sum::<f64>();
    }

    solve_linear_system(normal, rhs)
}

fn solve_linear_system(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    if matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
        return None;
    }

    for pivot_index in 0..n {
        let (best_row, best_value) = (pivot_index..n)
            .map(|row| (row, matrix[row][pivot_index].abs()))
            .max_by(|left, right| {
                left.1
                    .partial_cmp(&right.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })?;
        if best_value <= f64::EPSILON {
            return None;
        }
        if best_row != pivot_index {
            matrix.swap(best_row, pivot_index);
            rhs.swap(best_row, pivot_index);
        }

        let pivot = matrix[pivot_index][pivot_index];
        for col in pivot_index..n {
            matrix[pivot_index][col] /= pivot;
        }
        rhs[pivot_index] /= pivot;

        for row in 0..n {
            if row == pivot_index {
                continue;
            }
            let factor = matrix[row][pivot_index];
            if factor.abs() <= f64::EPSILON {
                continue;
            }
            for col in pivot_index..n {
                matrix[row][col] -= factor * matrix[pivot_index][col];
            }
            rhs[row] -= factor * rhs[pivot_index];
        }
    }

    Some(rhs)
}

fn monotone_cubic_tangents(x: &[f64], y: &[f64]) -> Option<Vec<f64>> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }

    let n = x.len();
    let mut h = Vec::with_capacity(n - 1);
    let mut delta = Vec::with_capacity(n - 1);
    for index in 0..n - 1 {
        let step = x[index + 1] - x[index];
        if step <= f64::EPSILON {
            return None;
        }
        h.push(step);
        delta.push((y[index + 1] - y[index]) / step);
    }

    let mut tangents = vec![0.0; n];
    tangents[0] = delta[0];
    tangents[n - 1] = delta[n - 2];

    for index in 1..n - 1 {
        if delta[index - 1] == 0.0
            || delta[index] == 0.0
            || delta[index - 1].signum() != delta[index].signum()
        {
            tangents[index] = 0.0;
            continue;
        }

        let w1 = 2.0 * h[index] + h[index - 1];
        let w2 = h[index] + 2.0 * h[index - 1];
        tangents[index] = (w1 + w2) / ((w1 / delta[index - 1]) + (w2 / delta[index]));
    }

    Some(tangents)
}

fn eval_monotone_cubic_spline(x: &[f64], y: &[f64], tangents: &[f64], xq: f64) -> f64 {
    if x.len() == 1 {
        return y[0];
    }
    if xq <= x[0] {
        return y[0] + tangents[0] * (xq - x[0]);
    }
    if xq >= x[x.len() - 1] {
        return y[y.len() - 1] + tangents[tangents.len() - 1] * (xq - x[x.len() - 1]);
    }

    let upper = x.partition_point(|value| *value <= xq);
    let index = upper.saturating_sub(1).min(x.len() - 2);
    let h = x[index + 1] - x[index];
    let t = (xq - x[index]) / h;
    let t2 = t * t;
    let t3 = t2 * t;

    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;

    h00 * y[index] + h10 * h * tangents[index] + h01 * y[index + 1] + h11 * h * tangents[index + 1]
}

fn eval_polynomial(coefficients: &[f64], x: f64) -> f64 {
    coefficients
        .iter()
        .enumerate()
        .map(|(power, coefficient)| coefficient * x.powi(power as i32))
        .sum::<f64>()
}

fn compute_ladder_qc_metrics(expected: &[f64], predicted: &[f64]) -> LadderQcMetrics {
    if expected.len() != predicted.len() || expected.is_empty() {
        return LadderQcMetrics {
            r2: f64::NEG_INFINITY,
            mean_abs_error_bp: f64::INFINITY,
            max_abs_error_bp: f64::INFINITY,
            monotonic_on_ladder: false,
        };
    }

    let mean_expected = expected.iter().sum::<f64>() / expected.len() as f64;
    let ss_tot = expected
        .iter()
        .map(|value| {
            let delta = value - mean_expected;
            delta * delta
        })
        .sum::<f64>();
    let residuals = expected
        .iter()
        .zip(predicted.iter())
        .map(|(exp, pred)| exp - pred)
        .collect::<Vec<_>>();
    let ss_res = residuals.iter().map(|value| value * value).sum::<f64>();
    let abs_errors = residuals
        .iter()
        .map(|value| value.abs())
        .collect::<Vec<_>>();
    let monotonic_on_ladder = predicted.windows(2).all(|window| window[1] > window[0]);

    LadderQcMetrics {
        r2: if ss_tot <= f64::EPSILON {
            f64::NEG_INFINITY
        } else {
            1.0 - (ss_res / ss_tot)
        },
        mean_abs_error_bp: abs_errors.iter().sum::<f64>() / abs_errors.len() as f64,
        max_abs_error_bp: abs_errors.into_iter().fold(0.0, f64::max),
        monotonic_on_ladder,
    }
}

fn build_sample_mapping_preview(
    sample_trace: &[f64],
    predicted_basepairs: &[f64],
) -> Option<SampleMappingPreview> {
    if sample_trace.is_empty() {
        return None;
    }
    if sample_trace.len() != predicted_basepairs.len() {
        return None;
    }

    let mut mapped = Vec::with_capacity(sample_trace.len());
    for (time, intensity) in sample_trace.iter().enumerate() {
        let basepair = predicted_basepairs[time];
        if !basepair.is_finite() || basepair < 0.0 {
            continue;
        }
        mapped.push(SampleMappingPoint {
            time,
            intensity: *intensity,
            basepair: round2(basepair),
        });
    }

    if mapped.is_empty() {
        return None;
    }

    let monotonic_unique = mapped
        .windows(2)
        .all(|window| window[1].basepair > window[0].basepair);
    let preview = sample_mapping_preview_points(&mapped);
    let sample_peak_preview = build_sample_peak_preview(sample_trace, &mapped);
    let assay_group_preview = build_assay_group_preview(&sample_peak_preview);
    let min_basepair = mapped.first().map(|point| point.basepair).unwrap_or(0.0);
    let max_basepair = mapped.last().map(|point| point.basepair).unwrap_or(0.0);

    Some(SampleMappingPreview {
        points_retained: mapped.len(),
        min_basepair,
        max_basepair,
        monotonic_unique,
        preview,
        sample_peak_preview,
        assay_group_preview,
    })
}

fn sample_mapping_preview_points(mapped: &[SampleMappingPoint]) -> Vec<SampleMappingPoint> {
    if mapped.len() <= 12 {
        return mapped.to_vec();
    }

    let mut preview = Vec::with_capacity(12);
    preview.extend_from_slice(&mapped[..4]);
    let mid = mapped.len() / 2;
    let middle_start = mid.saturating_sub(2);
    let middle_end = (middle_start + 4).min(mapped.len());
    preview.extend_from_slice(&mapped[middle_start..middle_end]);
    preview.extend_from_slice(&mapped[mapped.len().saturating_sub(4)..]);
    preview
}

fn round2(value: f64) -> f64 {
    (value * 100.0).round() / 100.0
}

fn build_sample_peak_preview(
    sample_trace: &[f64],
    mapped: &[SampleMappingPoint],
) -> Vec<SamplePeakPreview> {
    if sample_trace.is_empty() || mapped.is_empty() {
        return Vec::new();
    }

    let min_height = estimate_sample_peak_min_height(sample_trace);
    let peaks = find_peaks(sample_trace, min_height, 8);
    let mut preview = peaks
        .into_iter()
        .filter_map(|peak| {
            mapped
                .binary_search_by_key(&peak.index, |point| point.time)
                .ok()
                .and_then(|mapped_index| mapped.get(mapped_index))
                .map(|point| SamplePeakPreview {
                    time: peak.index,
                    intensity: round2(peak.height),
                    basepair: point.basepair,
                })
        })
        .collect::<Vec<_>>();

    preview.sort_by(|left, right| {
        right
            .intensity
            .partial_cmp(&left.intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                left.basepair
                    .partial_cmp(&right.basepair)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    if preview.len() > 12 {
        preview.truncate(12);
    }
    preview.sort_by(|left, right| {
        left.basepair
            .partial_cmp(&right.basepair)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    preview
}

fn estimate_sample_peak_min_height(sample_trace: &[f64]) -> f64 {
    let mut positives = sample_trace
        .iter()
        .copied()
        .filter(|value| *value > 0.0)
        .collect::<Vec<_>>();
    if positives.is_empty() {
        return 1.0;
    }
    positives.sort_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
    let median = positives[positives.len() / 2];
    let max_value = positives.last().copied().unwrap_or(median);
    let adaptive_floor = (max_value * 0.05).max(median * 3.0);
    adaptive_floor.max(25.0)
}

fn build_assay_group_preview(peaks: &[SamplePeakPreview]) -> Vec<SamplePeakGroupPreview> {
    if peaks.is_empty() {
        return Vec::new();
    }

    let mut sorted = peaks.to_vec();
    sorted.sort_by(|left, right| {
        left.basepair
            .partial_cmp(&right.basepair)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut groups: Vec<Vec<SamplePeakPreview>> = Vec::new();
    for peak in sorted {
        let should_start_new = groups
            .last()
            .and_then(|group| group.last())
            .map(|last| peak.basepair - last.basepair > SAMPLE_ASSAY_GROUP_DISTANCE_BP)
            .unwrap_or(true);
        if should_start_new {
            groups.push(vec![peak]);
        } else if let Some(group) = groups.last_mut() {
            group.push(peak);
        }
    }

    groups
        .into_iter()
        .enumerate()
        .filter_map(|(index, group)| {
            let max_intensity = group
                .iter()
                .map(|peak| peak.intensity)
                .fold(f64::NEG_INFINITY, f64::max);
            if !max_intensity.is_finite() || max_intensity <= 0.0 {
                return None;
            }

            let kept = group
                .into_iter()
                .filter(|peak| peak.intensity / max_intensity >= SAMPLE_ASSAY_MIN_RATIO)
                .collect::<Vec<_>>();
            if kept.is_empty() {
                return None;
            }

            let start_basepair = kept.first().map(|peak| peak.basepair).unwrap_or(0.0);
            let end_basepair = kept
                .last()
                .map(|peak| peak.basepair)
                .unwrap_or(start_basepair);
            let cluster_width_bp = round2((end_basepair - start_basepair).max(0.0));
            let dominant_peak = kept
                .iter()
                .max_by(|left, right| {
                    left.intensity
                        .partial_cmp(&right.intensity)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .cloned()?;
            let mut sorted_by_intensity = kept.clone();
            sorted_by_intensity.sort_by(|left, right| {
                right
                    .intensity
                    .partial_cmp(&left.intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let dominant_ratio_vs_second = sorted_by_intensity
                .get(1)
                .map(|second| dominant_peak.intensity / second.intensity.max(1.0));
            let clonal_candidate = kept.len() <= CLONAL_MAX_LABELLED_PEAKS
                && dominant_ratio_vs_second
                    .map(|ratio| ratio >= CLONAL_DOMINANCE_RATIO)
                    .unwrap_or(true);
            Some(SamplePeakGroupPreview {
                group_id: index + 1,
                start_basepair,
                end_basepair,
                cluster_width_bp,
                max_intensity: round2(max_intensity),
                dominant_peak_basepair: round2(dominant_peak.basepair),
                dominant_peak_intensity: round2(dominant_peak.intensity),
                dominant_ratio_vs_second: dominant_ratio_vs_second.map(round2),
                kept_peak_count: kept.len(),
                clonal_candidate,
                peaks: kept,
            })
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct ClonalityAssayDef {
    name: &'static str,
    channels: &'static [&'static str],
    bp_min: f64,
    bp_max: f64,
    aliases: &'static [&'static str],
}

const CLONALITY_ASSAYS: &[ClonalityAssayDef] = &[
    ClonalityAssayDef {
        name: "FR1",
        channels: &["DATA1"],
        bp_min: 280.0,
        bp_max: 420.0,
        aliases: &["FR1"],
    },
    ClonalityAssayDef {
        name: "FR2",
        channels: &["DATA1"],
        bp_min: 200.0,
        bp_max: 400.0,
        aliases: &["FR2"],
    },
    ClonalityAssayDef {
        name: "FR3",
        channels: &["DATA2"],
        bp_min: 60.0,
        bp_max: 220.0,
        aliases: &["FR3"],
    },
    ClonalityAssayDef {
        name: "IGK",
        channels: &["DATA1", "DATA2"],
        bp_min: 90.0,
        bp_max: 330.0,
        aliases: &["IGK"],
    },
    ClonalityAssayDef {
        name: "KDE",
        channels: &["DATA3"],
        bp_min: 190.0,
        bp_max: 410.0,
        aliases: &["KDE"],
    },
    ClonalityAssayDef {
        name: "TCRgA",
        channels: &["DATA1", "DATA2"],
        bp_min: 110.0,
        bp_max: 290.0,
        aliases: &["TCRGA", "TCRG"],
    },
    ClonalityAssayDef {
        name: "TCRgB",
        channels: &["DATA1", "DATA2"],
        bp_min: 60.0,
        bp_max: 250.0,
        aliases: &["TCRGB", "TCRG"],
    },
    ClonalityAssayDef {
        name: "DHJH_D",
        channels: &["DATA2"],
        bp_min: 90.0,
        bp_max: 440.0,
        aliases: &["DHJHD", "DHJH"],
    },
    ClonalityAssayDef {
        name: "DHJH_E",
        channels: &["DATA1"],
        bp_min: 65.0,
        bp_max: 160.0,
        aliases: &["DHJHE", "DHJH"],
    },
    ClonalityAssayDef {
        name: "TCRbA",
        channels: &["DATA1", "DATA2"],
        bp_min: 210.0,
        bp_max: 310.0,
        aliases: &["TCRBA", "TCRB"],
    },
    ClonalityAssayDef {
        name: "TCRbB",
        channels: &["DATA1", "DATA2"],
        bp_min: 210.0,
        bp_max: 310.0,
        aliases: &["TCRBB", "TCRB"],
    },
    ClonalityAssayDef {
        name: "TCRbC",
        channels: &["DATA1", "DATA2"],
        bp_min: 140.0,
        bp_max: 360.0,
        aliases: &["TCRBC", "TCRB"],
    },
];

const FLT3_ASSAYS: &[Flt3AssayDef] = &[
    Flt3AssayDef {
        name: "FLT3-ITD",
        channels: &["DATA1", "DATA2"],
        bp_min: 50.0,
        bp_max: 1000.0,
        wt_bp: 330.0,
        wt_tolerance_bp: 8.0,
        mutant_bp_min: 335.0,
        mutant_bp_max: 1000.0,
        positive_ratio: 0.02,
        aliases: &["FLT3ITD", "ITD", "ITDR", "RATIO"],
    },
    Flt3AssayDef {
        name: "FLT3-D835",
        channels: &["DATA3"],
        bp_min: 50.0,
        bp_max: 250.0,
        wt_bp: 80.0,
        wt_tolerance_bp: 4.0,
        mutant_bp_min: 121.0,
        mutant_bp_max: 130.5,
        positive_ratio: 0.05,
        aliases: &["FLT3D835", "D835", "TKD", "KUTTING"],
    },
    Flt3AssayDef {
        name: "NPM1",
        channels: &["DATA3"],
        bp_min: 50.0,
        bp_max: 1000.0,
        wt_bp: 300.0,
        wt_tolerance_bp: 3.0,
        mutant_bp_min: 303.0,
        mutant_bp_max: 305.0,
        positive_ratio: 0.01,
        aliases: &["NPM1", "NPM"],
    },
];

fn build_clonality_preview(
    file_name: &str,
    sample_channel: &str,
    assay_groups: &[SamplePeakGroupPreview],
) -> ClonalityPreview {
    let normalized_file_name = normalize_assay_token(file_name);
    let mut ranked_assays = CLONALITY_ASSAYS
        .iter()
        .filter_map(|assay| {
            let compatible_channel = assay.channels.contains(&sample_channel);
            let matched_by_filename = assay
                .aliases
                .iter()
                .any(|alias| normalized_file_name.contains(&normalize_assay_token(alias)));

            let matched_groups = assay_groups
                .iter()
                .filter_map(|group| {
                    let overlap_start = group.start_basepair.max(assay.bp_min);
                    let overlap_end = group.end_basepair.min(assay.bp_max);
                    if overlap_end < overlap_start {
                        return None;
                    }
                    Some(ClonalityGroupMatch {
                        group_id: group.group_id,
                        overlap_start_bp: round2(overlap_start),
                        overlap_end_bp: round2(overlap_end),
                        peak_count: group.kept_peak_count,
                        cluster_width_bp: group.cluster_width_bp,
                        dominant_peak_basepair: group.dominant_peak_basepair,
                        dominant_peak_intensity: group.dominant_peak_intensity,
                        dominant_ratio_vs_second: group.dominant_ratio_vs_second,
                        clonal_candidate: group.clonal_candidate,
                    })
                })
                .collect::<Vec<_>>();

            if !matched_by_filename && matched_groups.is_empty() && !compatible_channel {
                return None;
            }

            let mut score = 0.0;
            if compatible_channel {
                score += 1.0;
            }
            if matched_by_filename {
                score += 3.0;
            }
            score += matched_groups.len() as f64 * 2.0;
            let clonal_group_count = matched_groups
                .iter()
                .filter(|group| group.clonal_candidate)
                .count();
            score += clonal_group_count as f64 * 1.5;
            let best_dominant_ratio = matched_groups
                .iter()
                .filter_map(|group| group.dominant_ratio_vs_second)
                .max_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
            if let Some(ratio) = best_dominant_ratio {
                score += (ratio - 1.0).max(0.0).min(4.0);
            }
            score += matched_groups
                .iter()
                .map(|group| group.overlap_end_bp - group.overlap_start_bp)
                .sum::<f64>()
                / 100.0;

            Some(ClonalityAssayMatch {
                assay_name: assay.name.to_owned(),
                matched_by_filename,
                compatible_channel,
                score: round2(score),
                clonal_group_count,
                best_dominant_ratio: best_dominant_ratio.map(round2),
                matched_groups,
            })
        })
        .collect::<Vec<_>>();

    ranked_assays.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.assay_name.cmp(&right.assay_name))
    });
    ranked_assays.truncate(5);

    ClonalityPreview {
        sample_channel: sample_channel.to_owned(),
        ranked_assays,
    }
}

fn detect_flt3_assay(file_name: &str) -> Option<&'static Flt3AssayDef> {
    let token = normalize_assay_token(file_name);
    FLT3_ASSAYS.iter().find(|assay| {
        assay
            .aliases
            .iter()
            .any(|alias| token.contains(&normalize_assay_token(alias)))
    })
}

fn build_flt3_preview(
    file_name: &str,
    sample_channel: &str,
    sample_peaks: &[SamplePeakPreview],
) -> Flt3Preview {
    let assay = detect_flt3_assay(file_name)
        .copied()
        .unwrap_or(FLT3_ASSAYS[0]);
    let compatible_channel = assay.channels.contains(&sample_channel);
    let matched_by_filename = detect_flt3_assay(file_name)
        .map(|detected| detected.name == assay.name)
        .unwrap_or(false);

    let assay_peaks = sample_peaks
        .iter()
        .filter(|peak| peak.basepair >= assay.bp_min && peak.basepair <= assay.bp_max)
        .cloned()
        .collect::<Vec<_>>();

    let wt_peak = assay_peaks
        .iter()
        .filter(|peak| (peak.basepair - assay.wt_bp).abs() <= assay.wt_tolerance_bp)
        .min_by(|left, right| {
            (left.basepair - assay.wt_bp)
                .abs()
                .partial_cmp(&(right.basepair - assay.wt_bp).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    right
                        .intensity
                        .partial_cmp(&left.intensity)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        })
        .cloned();

    let mut mutant_peaks = assay_peaks
        .into_iter()
        .filter(|peak| peak.basepair >= assay.mutant_bp_min && peak.basepair <= assay.mutant_bp_max)
        .collect::<Vec<_>>();
    mutant_peaks.sort_by(|left, right| {
        right
            .intensity
            .partial_cmp(&left.intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if mutant_peaks.len() > 3 {
        mutant_peaks.truncate(3);
    }

    let strongest_mutant_ratio = wt_peak.as_ref().and_then(|wt| {
        mutant_peaks
            .first()
            .map(|mutant| round2(mutant.intensity / wt.intensity.max(1.0)))
    });
    let positive_call = if let Some(ratio) = strongest_mutant_ratio {
        ratio >= assay.positive_ratio && !mutant_peaks.is_empty()
    } else {
        !mutant_peaks.is_empty()
    };

    Flt3Preview {
        assay_name: assay.name.to_owned(),
        matched_by_filename,
        compatible_channel,
        assay_bp_min: assay.bp_min,
        assay_bp_max: assay.bp_max,
        wt_peak,
        mutant_peaks,
        strongest_mutant_ratio,
        positive_call,
    }
}

fn normalize_assay_token(value: &str) -> String {
    value
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .flat_map(|ch| {
            ch.to_ascii_uppercase()
                .to_string()
                .chars()
                .collect::<Vec<_>>()
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
struct RefinementCandidate {
    changed_step_indices: Vec<usize>,
    original_scan_indices: Vec<usize>,
    refined_scan_indices: Vec<usize>,
    sizing_model: SizingModelPreview,
}

fn refine_best_combination(
    peak_pool: &[usize],
    current_scans: &[usize],
    ladder_sizes: &[f64],
    current_model: &SizingModelPreview,
) -> Option<RefinementCandidate> {
    if current_scans.len() != ladder_sizes.len() || current_scans.len() < 3 {
        return None;
    }

    let residuals = ladder_sizes
        .iter()
        .zip(current_model.predicted_ladder_basepairs.iter())
        .map(|(expected, predicted)| expected - predicted)
        .collect::<Vec<_>>();
    let mut ranked_steps = residuals
        .iter()
        .enumerate()
        .map(|(index, residual)| (index, residual.abs()))
        .filter(|(_, abs_residual)| *abs_residual >= MIN_REFINEMENT_TRIGGER_BP)
        .collect::<Vec<_>>();
    ranked_steps.sort_by(|left, right| {
        right
            .1
            .partial_cmp(&left.1)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    ranked_steps.truncate(MAX_REFINEMENT_STEPS);
    if ranked_steps.is_empty() {
        return None;
    }

    let mut option_buckets = Vec::new();
    let mut changed_indices = Vec::new();

    for (step_index, _) in ranked_steps {
        let current_scan = current_scans[step_index] as f64;
        let slope = local_bp_per_scan(
            current_scans,
            &current_model.predicted_ladder_basepairs,
            step_index,
        );
        if !slope.is_finite() || slope.abs() <= 1e-6 {
            continue;
        }
        let target_scan = current_scan + (residuals[step_index] / slope);
        let lower_bound = if step_index == 0 {
            f64::NEG_INFINITY
        } else {
            current_scans[step_index - 1] as f64 + 6.0
        };
        let upper_bound = if step_index + 1 >= current_scans.len() {
            f64::INFINITY
        } else {
            current_scans[step_index + 1] as f64 - 6.0
        };
        let options = refinement_options(
            peak_pool,
            current_scan,
            target_scan,
            lower_bound,
            upper_bound,
        );
        if options.len() < 2 {
            continue;
        }
        changed_indices.push(step_index);
        option_buckets.push(options);
    }

    if option_buckets.is_empty() {
        return None;
    }

    let baseline_score = model_score(&current_model.qc_metrics);
    let mut best_candidate: Option<RefinementCandidate> = None;
    let mut best_score = baseline_score;
    let mut current_trial = current_scans.to_vec();
    try_refinement_combinations(
        0,
        &changed_indices,
        &option_buckets,
        &mut current_trial,
        current_scans,
        ladder_sizes,
        &mut best_candidate,
        &mut best_score,
    );
    best_candidate
}

fn try_refinement_combinations(
    bucket_index: usize,
    changed_indices: &[usize],
    option_buckets: &[Vec<usize>],
    trial_scans: &mut Vec<usize>,
    original_scans: &[usize],
    ladder_sizes: &[f64],
    best_candidate: &mut Option<RefinementCandidate>,
    best_score: &mut (f64, f64, f64),
) {
    if bucket_index == changed_indices.len() {
        if trial_scans == original_scans {
            return;
        }
        if !trial_scans.windows(2).all(|window| window[1] > window[0]) {
            return;
        }
        let Some(model) = fit_best_sizing_model(trial_scans, ladder_sizes, &[]) else {
            return;
        };
        let score = model_score(&model.qc_metrics);
        if score < *best_score {
            *best_score = score;
            *best_candidate = Some(RefinementCandidate {
                changed_step_indices: changed_indices.to_vec(),
                original_scan_indices: original_scans.to_vec(),
                refined_scan_indices: trial_scans.clone(),
                sizing_model: model,
            });
        }
        return;
    }

    let step_index = changed_indices[bucket_index];
    let original_value = trial_scans[step_index];
    for candidate in &option_buckets[bucket_index] {
        trial_scans[step_index] = *candidate;
        if step_index > 0 && trial_scans[step_index] <= trial_scans[step_index - 1] {
            continue;
        }
        if step_index + 1 < trial_scans.len()
            && trial_scans[step_index] >= trial_scans[step_index + 1]
        {
            continue;
        }
        try_refinement_combinations(
            bucket_index + 1,
            changed_indices,
            option_buckets,
            trial_scans,
            original_scans,
            ladder_sizes,
            best_candidate,
            best_score,
        );
    }
    trial_scans[step_index] = original_value;
}

fn refinement_options(
    peak_pool: &[usize],
    current_scan: f64,
    target_scan: f64,
    lower_bound: f64,
    upper_bound: f64,
) -> Vec<usize> {
    let mut options = peak_pool
        .iter()
        .copied()
        .filter(|scan| {
            let scan_f64 = *scan as f64;
            scan_f64 >= lower_bound
                && scan_f64 <= upper_bound
                && (scan_f64 - current_scan).abs() <= MAX_REFINEMENT_RADIUS_SCANS
                || (scan_f64 - target_scan).abs() <= MAX_REFINEMENT_RADIUS_SCANS
        })
        .collect::<Vec<_>>();

    options.sort_by(|left, right| {
        let left_score = (
            (*left as f64 - target_scan).abs(),
            (*left as f64 - current_scan).abs(),
        );
        let right_score = (
            (*right as f64 - target_scan).abs(),
            (*right as f64 - current_scan).abs(),
        );
        left_score
            .partial_cmp(&right_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    options.dedup();

    let current_scan_usize = current_scan.round() as usize;
    if !options.contains(&current_scan_usize) {
        options.insert(0, current_scan_usize);
    }
    options.truncate(MAX_REFINEMENT_OPTIONS_PER_STEP);
    if options.is_empty() {
        vec![current_scan_usize]
    } else {
        options
    }
}

fn local_bp_per_scan(scans: &[usize], predicted_bps: &[f64], index: usize) -> f64 {
    if scans.len() != predicted_bps.len() || scans.len() < 2 {
        return f64::NAN;
    }

    if index == 0 {
        let dx = scans[1] as f64 - scans[0] as f64;
        return (predicted_bps[1] - predicted_bps[0]) / dx;
    }
    if index + 1 >= scans.len() {
        let last = scans.len() - 1;
        let dx = scans[last] as f64 - scans[last - 1] as f64;
        return (predicted_bps[last] - predicted_bps[last - 1]) / dx;
    }

    let dx = scans[index + 1] as f64 - scans[index - 1] as f64;
    (predicted_bps[index + 1] - predicted_bps[index - 1]) / dx
}

fn model_score(metrics: &LadderQcMetrics) -> (f64, f64, f64) {
    (
        -metrics.r2,
        metrics.max_abs_error_bp,
        metrics.mean_abs_error_bp,
    )
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::ladders::LadderKind;
    use crate::signal::Peak;

    use super::{
        SamplePeakGroupPreview, SamplePeakPreview, SizingModelPreview, build_assay_group_preview,
        build_clonality_preview, build_flt3_preview, build_sample_mapping_preview,
        compute_ladder_qc_metrics, curvature_score, estimate_combination_count_capped,
        eval_polynomial, fit_polynomial_least_squares, generate_peak_combinations,
        ladder_peak_sequence_penalty, quadratic_fit_r2, refine_best_combination,
        select_best_combination, select_ladder_peaks,
    };

    #[test]
    fn quadratic_fit_scores_perfect_curve() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![1.0, 3.0, 7.0, 13.0, 21.0];
        let r2 = quadratic_fit_r2(&x, &y);
        assert!(r2 > 0.999);
    }

    #[test]
    fn generate_peak_combinations_respects_gap_and_target_length() {
        let peaks = vec![100, 205, 310, 415, 900];
        let combinations = generate_peak_combinations(&peaks, 4, 130, 100);
        assert_eq!(combinations, vec![vec![100, 205, 310, 415]]);
    }

    #[test]
    fn combination_estimator_caps_large_search_spaces() {
        let peaks = (0..30).map(|index| 100 + index * 10).collect::<Vec<_>>();
        let count = estimate_combination_count_capped(&peaks, 10, 40, 500);
        assert_eq!(count, 500);
    }

    #[test]
    fn select_best_combination_picks_lowest_curvature_candidate() {
        let ladder_sizes = vec![35.0, 50.0, 75.0, 100.0];
        let combinations = vec![vec![100, 200, 300, 400], vec![100, 200, 380, 620]];
        let best = select_best_combination(
            &combinations,
            &ladder_sizes,
            LadderKind::Liz500250,
            &[],
        )
            .expect("a best combination should be selected");
        assert_eq!(best.indices, vec![100, 200, 380, 620]);
        assert!(curvature_score(&ladder_sizes, &best.indices).is_finite());
    }

    #[test]
    fn select_ladder_peaks_prefers_raw_trace_before_baseline_fallback() {
        let raw = vec![0.0, 10.0, 400.0, 20.0, 0.0, 10.0, 450.0, 10.0, 0.0];
        let corrected = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let peaks = select_ladder_peaks(&raw, &corrected, &corrected, 100.0, 60.0, 2, 10, 4);
        let indices = peaks.iter().map(|peak| peak.index).collect::<Vec<_>>();
        assert_eq!(indices, vec![2, 6]);
    }

    #[test]
    fn select_ladder_peaks_prefers_corrected_when_raw_has_only_blob_peaks() {
        let raw = vec![0.0, 0.0, 900.0, 0.0, 0.0, 850.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let corrected = vec![
            0.0, 120.0, 0.0, 130.0, 0.0, 125.0, 0.0, 140.0, 0.0, 135.0, 0.0,
        ];
        let peaks =
            select_ladder_peaks(&raw, &corrected, &corrected, 100.0, 60.0, 1, 10, 5);
        let indices = peaks.iter().map(|peak| peak.index).collect::<Vec<_>>();
        for expected in [1usize, 3, 5, 7, 9] {
            assert!(indices.contains(&expected));
        }
    }

    #[test]
    fn ladder_domain_penalty_applies_early_anchor_checks_for_liz() {
        let peak = Peak {
            index: 2100,
            height: 120.0,
            prominence: 100.0,
            width: 5.0,
            local_baseline: 0.0,
            score: 120.0,
        };
        let peak_map = [(2100usize, peak.clone())].into_iter().collect();
        let penalty = super::ladder_domain_penalty(
            LadderKind::Liz500250,
            &[2100, 2300, 2500, 2700],
            &peak_map,
            &[peak],
        );
        assert!(penalty > 0.0);
    }

    #[test]
    fn top_peak_candidates_prefers_scored_peak_shape_over_raw_height() {
        let values = vec![0.0, 160.0, 0.0, 70.0, 110.0, 130.0, 110.0, 70.0, 0.0];
        let peaks = super::top_peak_candidates(&values, 50.0, 1, 1);
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].index, 5);
    }

    #[test]
    fn ladder_peak_sequence_penalty_dislikes_crowded_dirty_start() {
        let ladder_sizes = LadderKind::Rox400Hd.sizes().to_vec();
        let clean = vec![
            1680usize, 1734, 1901, 1960, 2073, 2250, 2310, 2428, 2488, 2548, 2671, 2791, 2915,
            3041, 3103, 3165, 3290, 3414, 3538, 3662, 3784,
        ];
        let crowded = vec![
            1516usize, 1640, 1901, 1960, 2073, 2250, 2310, 2428, 2488, 2548, 2671, 2791, 2915,
            3041, 3103, 3165, 3290, 3414, 3538, 3662, 3784,
        ];
        let mut peak_map = BTreeMap::new();
        for scan in &clean {
            peak_map.insert(
                *scan,
                Peak {
                    index: *scan,
                    height: 240.0,
                    prominence: 220.0,
                    width: 4.0,
                    local_baseline: 10.0,
                    score: 260.0,
                },
            );
        }
        peak_map.insert(
            1516,
            Peak {
                index: 1516,
                height: 260.0,
                prominence: 80.0,
                width: 12.0,
                local_baseline: 150.0,
                score: 110.0,
            },
        );
        peak_map.insert(
            1640,
            Peak {
                index: 1640,
                height: 230.0,
                prominence: 120.0,
                width: 8.0,
                local_baseline: 90.0,
                score: 135.0,
            },
        );
        let peak_features = peak_map.values().cloned().collect::<Vec<_>>();
        let clean_penalty =
            ladder_peak_sequence_penalty(&clean, &ladder_sizes, &peak_map, &peak_features);
        let crowded_penalty =
            ladder_peak_sequence_penalty(&crowded, &ladder_sizes, &peak_map, &peak_features);
        assert!(crowded_penalty > clean_penalty);
    }

    #[test]
    fn polynomial_fit_recovers_linear_mapping() {
        let x = vec![100.0, 200.0, 300.0, 400.0];
        let y = vec![50.0, 100.0, 150.0, 200.0];
        let coeffs =
            fit_polynomial_least_squares(&x, &y, 1).expect("linear least squares should succeed");
        let predicted = x
            .iter()
            .map(|value| eval_polynomial(&coeffs, *value))
            .collect::<Vec<_>>();
        let qc = compute_ladder_qc_metrics(&y, &predicted);
        assert!(qc.r2 > 0.999999);
        assert!(qc.mean_abs_error_bp < 1e-6);
        assert!(qc.monotonic_on_ladder);
    }

    #[test]
    fn local_refinement_can_improve_a_nearby_candidate() {
        let ladder_sizes = vec![35.0, 50.0, 75.0, 100.0];
        let scans = vec![100, 200, 320, 405];
        let model = SizingModelPreview {
            strategy: "test".to_owned(),
            degree: 1,
            coefficients: fit_polynomial_least_squares(
                &scans.iter().map(|value| *value as f64).collect::<Vec<_>>(),
                &ladder_sizes,
                1,
            )
            .expect("fit should work"),
            predicted_ladder_basepairs: scans.iter().map(|value| *value as f64).collect::<Vec<_>>(),
            qc_metrics: compute_ladder_qc_metrics(&ladder_sizes, &[35.0, 50.0, 80.0, 101.5]),
            sample_mapping: None,
        };
        let peak_pool = vec![100, 200, 300, 320, 405];
        let refined = refine_best_combination(&peak_pool, &scans, &ladder_sizes, &model)
            .expect("refinement should find an improvement");
        assert_eq!(refined.refined_scan_indices, vec![100, 200, 300, 405]);
        assert!(
            refined.sizing_model.qc_metrics.max_abs_error_bp < model.qc_metrics.max_abs_error_bp
        );
    }

    #[test]
    fn sample_mapping_preview_filters_negative_basepairs_and_stays_monotonic() {
        let trace = vec![10.0, 20.0, 30.0, 40.0, 50.0];
        let predicted = vec![-1.0, 1.0, 3.0, 5.0, 7.0];
        let preview = build_sample_mapping_preview(&trace, &predicted)
            .expect("sample mapping preview should be created");
        assert_eq!(preview.points_retained, 4);
        assert_eq!(preview.min_basepair, 1.0);
        assert!(preview.monotonic_unique);
    }

    #[test]
    fn sample_mapping_preview_includes_detected_sample_peaks() {
        let trace = vec![
            0.0, 10.0, 200.0, 15.0, 0.0, 5.0, 0.0, 4.0, 0.0, 20.0, 250.0, 10.0, 0.0,
        ];
        let predicted = (0..trace.len())
            .map(|value| value as f64)
            .collect::<Vec<_>>();
        let preview = build_sample_mapping_preview(&trace, &predicted)
            .expect("sample mapping preview should be created");
        assert_eq!(preview.sample_peak_preview.len(), 2);
        assert_eq!(preview.sample_peak_preview[0].time, 2);
        assert_eq!(preview.sample_peak_preview[1].time, 10);
    }

    #[test]
    fn assay_group_preview_clusters_nearby_peaks_and_filters_weak_ones() {
        let peaks = vec![
            SamplePeakPreview {
                time: 10,
                intensity: 100.0,
                basepair: 100.0,
            },
            SamplePeakPreview {
                time: 12,
                intensity: 20.0,
                basepair: 104.0,
            },
            SamplePeakPreview {
                time: 20,
                intensity: 50.0,
                basepair: 130.0,
            },
            SamplePeakPreview {
                time: 22,
                intensity: 30.0,
                basepair: 138.0,
            },
        ];
        let groups = build_assay_group_preview(&peaks);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].kept_peak_count, 1);
        assert_eq!(groups[1].kept_peak_count, 2);
        assert!(groups[0].clonal_candidate);
        assert_eq!(groups[0].dominant_ratio_vs_second, None);
        assert_eq!(groups[1].dominant_ratio_vs_second, Some(1.67));
        assert!(!groups[1].clonal_candidate);
    }

    #[test]
    fn clonality_preview_prefers_filename_and_group_overlap_matches() {
        let groups = vec![SamplePeakGroupPreview {
            group_id: 1,
            start_basepair: 120.0,
            end_basepair: 150.0,
            cluster_width_bp: 30.0,
            max_intensity: 3000.0,
            dominant_peak_basepair: 140.0,
            dominant_peak_intensity: 3000.0,
            dominant_ratio_vs_second: Some(2.2),
            kept_peak_count: 1,
            clonal_candidate: true,
            peaks: vec![SamplePeakPreview {
                time: 100,
                intensity: 3000.0,
                basepair: 140.0,
            }],
        }];
        let preview =
            build_clonality_preview("26OUM04817_IGK_270326_B05_H9H1DI2F.fsa", "DATA1", &groups);
        assert!(!preview.ranked_assays.is_empty());
        assert_eq!(preview.ranked_assays[0].assay_name, "IGK");
        assert!(preview.ranked_assays[0].matched_by_filename);
        assert_eq!(preview.ranked_assays[0].clonal_group_count, 1);
        assert_eq!(preview.ranked_assays[0].best_dominant_ratio, Some(2.2));
        assert!(preview.ranked_assays[0].matched_groups[0].clonal_candidate);
    }

    #[test]
    fn flt3_preview_detects_itd_mutant_peaks_from_filename() {
        let peaks = vec![
            SamplePeakPreview {
                time: 100,
                intensity: 10000.0,
                basepair: 330.2,
            },
            SamplePeakPreview {
                time: 120,
                intensity: 850.0,
                basepair: 351.0,
            },
            SamplePeakPreview {
                time: 140,
                intensity: 410.0,
                basepair: 371.2,
            },
        ];

        let preview = build_flt3_preview("ivs0000_flt3_itd_p1.fsa", "DATA1", &peaks);
        assert_eq!(preview.assay_name, "FLT3-ITD");
        assert!(preview.matched_by_filename);
        assert!(preview.compatible_channel);
        assert!(preview.positive_call);
        assert_eq!(
            preview.wt_peak.as_ref().map(|peak| peak.basepair),
            Some(330.2)
        );
        assert_eq!(preview.mutant_peaks.len(), 2);
        assert!(preview.strongest_mutant_ratio.is_some());
    }
}
