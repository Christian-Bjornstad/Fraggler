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
            let sample_mapping = build_sample_mapping_preview(sample_trace, &coefficients);
            return Some(SizingModelPreview {
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
    coefficients: &[f64],
) -> Option<SampleMappingPreview> {
    if sample_trace.is_empty() {
        return None;
    }

    let mut mapped = Vec::with_capacity(sample_trace.len());
    for (time, intensity) in sample_trace.iter().enumerate() {
        let basepair = eval_polynomial(coefficients, time as f64);
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
    let min_basepair = mapped.first().map(|point| point.basepair).unwrap_or(0.0);
    let max_basepair = mapped.last().map(|point| point.basepair).unwrap_or(0.0);

    Some(SampleMappingPreview {
        points_retained: mapped.len(),
        min_basepair,
        max_basepair,
        monotonic_unique,
        preview,
        sample_peak_preview,
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

    let derivative_coeffs = derivative_coefficients(&current_model.coefficients);
    let mut option_buckets = Vec::new();
    let mut changed_indices = Vec::new();

    for (step_index, _) in ranked_steps {
        let current_scan = current_scans[step_index] as f64;
        let derivative = eval_polynomial(&derivative_coeffs, current_scan);
        if !derivative.is_finite() || derivative.abs() <= 1e-6 {
            continue;
        }
        let target_scan = current_scan + (residuals[step_index] / derivative);
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

fn derivative_coefficients(coefficients: &[f64]) -> Vec<f64> {
    coefficients
        .iter()
        .enumerate()
        .skip(1)
        .map(|(power, coefficient)| *coefficient * power as f64)
        .collect()
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
        SizingModelPreview, build_sample_mapping_preview, compute_ladder_qc_metrics,
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
        let coeffs = vec![-1.0, 2.0];
        let preview = build_sample_mapping_preview(&trace, &coeffs)
            .expect("sample mapping preview should be created");
        assert_eq!(preview.points_retained, 4);
        assert_eq!(preview.min_basepair, 1.0);
        assert!(preview.monotonic_unique);
    }

    #[test]
    fn sample_mapping_preview_includes_detected_sample_peaks() {
        let trace = vec![0.0, 10.0, 200.0, 15.0, 0.0, 5.0, 0.0, 4.0, 0.0, 20.0, 250.0, 10.0, 0.0];
        let coeffs = vec![0.0, 1.0];
        let preview = build_sample_mapping_preview(&trace, &coeffs)
            .expect("sample mapping preview should be created");
        assert_eq!(preview.sample_peak_preview.len(), 2);
        assert_eq!(preview.sample_peak_preview[0].time, 2);
        assert_eq!(preview.sample_peak_preview[1].time, 10);
    }
}
