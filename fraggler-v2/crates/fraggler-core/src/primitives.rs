use std::collections::BTreeMap;

use camino::Utf8Path;
use serde::{Deserialize, Serialize};

use crate::abif::AbifRecord;
use crate::contract::AnalysisKind;
use crate::engine::EngineError;
use crate::ladders::LadderKind;
use crate::signal::{Peak, baseline_correct_nonnegative, find_peaks};

const MAX_CANDIDATE_COMBINATIONS: usize = 100_000;
const LADDER_MAX_GAP_EXPANSIONS: usize = 15;
const LADDER_GAP_EXPANSION_STEP: usize = 10;
const MAX_REFINEMENT_STEPS: usize = 3;
const MAX_REFINEMENT_OPTIONS_PER_STEP: usize = 5;
const MIN_REFINEMENT_TRIGGER_BP: f64 = 0.75;
const MAX_REFINEMENT_RADIUS_SCANS: f64 = 120.0;
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
    let sample_channel = data_channels
        .iter()
        .find(|channel| channel.as_str() != size_standard_channel)
        .cloned()
        .unwrap_or_else(|| size_standard_channel.clone());

    let ladder = suggested_ladder_kind(&record, &size_standard_channel, analysis_kind);
    let min_height = match ladder {
        LadderKind::Liz500250 => 300.0,
        LadderKind::Rox400Hd => 120.0,
        LadderKind::Gs500Rox => 120.0,
    };
    let min_distance = match ladder {
        LadderKind::Liz500250 => 15,
        LadderKind::Rox400Hd => 8,
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
    let corrected = baseline_correct_nonnegative(&size_standard_trace, 0.01, 100.0, 50)?;
    let ladder_peaks = select_ladder_peaks(
        &size_standard_trace,
        &corrected,
        min_height,
        min_distance,
        ladder.expected_peak_count() + 15,
    );
    let file_name = record.file_name.clone();
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

fn select_ladder_peaks(
    raw_trace: &[f64],
    corrected_trace: &[f64],
    min_height: f64,
    min_distance: usize,
    max_peaks: usize,
) -> Vec<Peak> {
    let raw_candidates = top_peak_candidates(raw_trace, min_height, min_distance, max_peaks);
    if !raw_candidates.is_empty() {
        return raw_candidates;
    }
    top_peak_candidates(corrected_trace, min_height, min_distance, max_peaks)
}

fn top_peak_candidates(
    values: &[f64],
    min_height: f64,
    min_distance: usize,
    max_peaks: usize,
) -> Vec<Peak> {
    let mut peaks = find_peaks(values, min_height, min_distance);
    peaks.sort_by(|left, right| {
        right
            .height
            .partial_cmp(&left.height)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.index.cmp(&right.index))
    });
    if peaks.len() > max_peaks {
        peaks.truncate(max_peaks);
    }
    peaks.sort_by_key(|peak| peak.index);
    peaks
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
    let mut max_allowed_peak_gap = estimate_max_allowed_peak_gap(&peak_indices, 2.0);
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

        combinations = generate_peak_combinations(
            &peak_indices,
            target_len,
            max_allowed_peak_gap,
            MAX_CANDIDATE_COMBINATIONS,
        );
        if !combinations.is_empty() {
            break;
        }
        max_allowed_peak_gap = max_allowed_peak_gap.saturating_add(LADDER_GAP_EXPANSION_STEP);
    }

    if candidate_generation_capped {
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
    let mut best = select_best_combination(&combinations, ladder_sizes);
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
            best = Some(CombinationScore {
                indices: refined.refined_scan_indices.clone(),
                curvature_score: refined_curvature,
                quadratic_r2: refined_qr2,
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

#[derive(Debug, Clone, PartialEq)]
struct CombinationScore {
    indices: Vec<usize>,
    curvature_score: f64,
    quadratic_r2: f64,
}

fn select_best_combination(
    combinations: &[Vec<usize>],
    ladder_sizes: &[f64],
) -> Option<CombinationScore> {
    combinations
        .iter()
        .filter(|combo| combo.len() == ladder_sizes.len())
        .map(|combo| CombinationScore {
            indices: combo.clone(),
            curvature_score: curvature_score(ladder_sizes, combo),
            quadratic_r2: quadratic_fit_r2(
                ladder_sizes,
                &combo.iter().map(|value| *value as f64).collect::<Vec<_>>(),
            ),
        })
        .min_by(|left, right| {
            left.curvature_score
                .partial_cmp(&right.curvature_score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    right
                        .quadratic_r2
                        .partial_cmp(&left.quadratic_r2)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| left.indices.cmp(&right.indices))
        })
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
    use super::{
        SamplePeakGroupPreview, SamplePeakPreview, SizingModelPreview, build_assay_group_preview,
        build_clonality_preview, build_sample_mapping_preview, compute_ladder_qc_metrics,
        curvature_score, estimate_combination_count_capped, eval_polynomial,
        fit_polynomial_least_squares, generate_peak_combinations, quadratic_fit_r2,
        refine_best_combination, select_best_combination, select_ladder_peaks,
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
        let best = select_best_combination(&combinations, &ladder_sizes)
            .expect("a best combination should be selected");
        assert_eq!(best.indices, vec![100, 200, 380, 620]);
        assert!(curvature_score(&ladder_sizes, &best.indices).is_finite());
    }

    #[test]
    fn select_ladder_peaks_prefers_raw_trace_before_baseline_fallback() {
        let raw = vec![0.0, 10.0, 400.0, 20.0, 0.0, 10.0, 450.0, 10.0, 0.0];
        let corrected = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let peaks = select_ladder_peaks(&raw, &corrected, 100.0, 2, 10);
        let indices = peaks.iter().map(|peak| peak.index).collect::<Vec<_>>();
        assert_eq!(indices, vec![2, 6]);
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
}
