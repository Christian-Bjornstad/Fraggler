use serde::{Deserialize, Serialize};

use crate::engine::EngineError;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Peak {
    pub index: usize,
    pub height: f64,
    pub prominence: f64,
    pub width: f64,
    pub local_baseline: f64,
    pub score: f64,
}

pub fn find_peaks(values: &[f64], min_height: f64, min_distance: usize) -> Vec<Peak> {
    if values.len() < 3 {
        return Vec::new();
    }

    let mut candidates = Vec::new();
    for index in 1..values.len() - 1 {
        let current = values[index];
        if current < min_height {
            continue;
        }
        let prev = values[index - 1];
        let next = values[index + 1];
        if current >= prev && current > next {
            let (local_baseline, prominence, width) = describe_peak_shape(values, index, current);
            candidates.push(Peak {
                index,
                height: current,
                prominence,
                width,
                local_baseline,
                score: peak_score(current, prominence, width, local_baseline),
            });
        }
    }

    candidates.sort_by(|left, right| {
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

    let mut accepted: Vec<Peak> = Vec::new();
    'candidate: for candidate in candidates {
        for kept in &accepted {
            let distance = candidate.index.abs_diff(kept.index);
            if distance < min_distance {
                continue 'candidate;
            }
        }
        accepted.push(candidate);
    }

    accepted.sort_by_key(|peak| peak.index);
    accepted
}

fn describe_peak_shape(values: &[f64], index: usize, current: f64) -> (f64, f64, f64) {
    let mut left_min = current;
    let mut left = index;
    while left > 0 {
        left -= 1;
        let value = values[left];
        if value < left_min {
            left_min = value;
        }
        if value > current {
            break;
        }
    }

    let mut right_min = current;
    let mut right = index;
    while right + 1 < values.len() {
        right += 1;
        let value = values[right];
        if value < right_min {
            right_min = value;
        }
        if value > current {
            break;
        }
    }

    let local_baseline = left_min.max(right_min);
    let prominence = (current - local_baseline).max(0.0);
    let half_height = local_baseline + 0.5 * prominence;

    let mut width_left = index;
    while width_left > 0 && values[width_left - 1] > half_height {
        width_left -= 1;
    }
    let mut width_right = index;
    while width_right + 1 < values.len() && values[width_right + 1] > half_height {
        width_right += 1;
    }
    let width = (width_right.saturating_sub(width_left) + 1) as f64;

    (local_baseline, prominence, width)
}

fn peak_score(height: f64, prominence: f64, width: f64, local_baseline: f64) -> f64 {
    let height_floor = height.max(1.0);
    let purity = (prominence / height_floor).clamp(0.0, 1.25);
    let width_term = if width <= 8.0 {
        width.clamp(1.0, 8.0).sqrt()
    } else {
        8.0_f64.sqrt() / (1.0 + 0.10 * (width - 8.0))
    };
    let baseline_ratio = (local_baseline.max(0.0) / height_floor).clamp(0.0, 1.5);
    let baseline_penalty = 1.0 / (1.0 + 1.6 * baseline_ratio);
    let purity_boost = 0.55 + 0.95 * purity;
    prominence * width_term * purity_boost * baseline_penalty + 0.08 * height
}

pub fn baseline_arpls(
    values: &[f64],
    ratio: f64,
    lam: f64,
    niter: usize,
) -> Result<Vec<f64>, EngineError> {
    if values.is_empty() {
        return Ok(Vec::new());
    }
    if values.len() < 3 {
        return Ok(values.to_vec());
    }

    let n = values.len();
    let mut weights = vec![1.0_f64; n];
    let mut baseline = vec![0.0_f64; n];
    let mut crit = f64::INFINITY;
    let target_ratio = if ratio <= 0.0 { 0.99 } else { ratio };
    let lambda = if lam <= 0.0 { 100.0 } else { lam };

    let (lower2_template, lower1_template, main_template, upper1_template, upper2_template) =
        difference_penalty_bands(n, lambda);

    let mut iterations = 0usize;
    while crit > target_ratio && iterations < niter {
        let main = main_template
            .iter()
            .zip(weights.iter())
            .map(|(penalty, weight)| penalty + weight)
            .collect::<Vec<_>>();
        let rhs = values
            .iter()
            .zip(weights.iter())
            .map(|(value, weight)| value * weight)
            .collect::<Vec<_>>();
        baseline = solve_pentadiagonal(
            &lower2_template,
            &lower1_template,
            &main,
            &upper1_template,
            &upper2_template,
            &rhs,
        )?;

        let residuals = values
            .iter()
            .zip(baseline.iter())
            .map(|(value, base)| value - base)
            .collect::<Vec<_>>();
        let negatives = residuals
            .iter()
            .copied()
            .filter(|value| *value < 0.0)
            .collect::<Vec<_>>();
        if negatives.is_empty() {
            break;
        }

        let mean = negatives.iter().sum::<f64>() / negatives.len() as f64;
        let variance = negatives
            .iter()
            .map(|value| {
                let delta = value - mean;
                delta * delta
            })
            .sum::<f64>()
            / negatives.len() as f64;
        let stddev = variance.sqrt();
        if !stddev.is_finite() || stddev <= f64::EPSILON {
            break;
        }

        let mut next_weights = Vec::with_capacity(n);
        for residual in &residuals {
            let exponent = 2.0 * (residual - (2.0 * stddev - mean)) / stddev;
            next_weights.push(1.0 / (1.0 + exponent.exp()));
        }

        let numerator = next_weights
            .iter()
            .zip(weights.iter())
            .map(|(next, current)| {
                let delta = next - current;
                delta * delta
            })
            .sum::<f64>()
            .sqrt();
        let denominator = weights
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        crit = if denominator > 0.0 {
            numerator / denominator
        } else {
            0.0
        };
        weights = next_weights;
        iterations += 1;
    }

    Ok(baseline)
}

pub fn baseline_correct_nonnegative(
    values: &[f64],
    ratio: f64,
    lam: f64,
    niter: usize,
) -> Result<Vec<f64>, EngineError> {
    let baseline = baseline_arpls(values, ratio, lam, niter)?;
    Ok(values
        .iter()
        .zip(baseline.iter())
        .map(|(value, base)| (value - base).max(0.0))
        .collect())
}

pub fn baseline_correct_guarded_nonnegative(
    values: &[f64],
    ratio: f64,
    lam: f64,
    niter: usize,
    bin_size: usize,
    quantile: f64,
) -> Result<Vec<f64>, EngineError> {
    if values.is_empty() {
        return Ok(Vec::new());
    }
    let envelope = rolling_quantile_baseline(values, bin_size, quantile);
    let baseline = baseline_arpls(values, ratio, lam, niter)?;
    let residual = values
        .iter()
        .zip(envelope.iter())
        .map(|(value, env)| value - env)
        .collect::<Vec<_>>();
    let residual_scale = stddev(&residual);
    let positive_excess = baseline
        .iter()
        .zip(envelope.iter())
        .map(|(base, env)| (base - env).max(0.0))
        .collect::<Vec<_>>();
    let mountain_score = slice_quantile(&positive_excess, 0.95);
    let upper_guard = 10.0_f64.max(0.10 * residual_scale);
    let lower_guard = 25.0_f64.max(2.5 * upper_guard);

    let guarded = if mountain_score <= upper_guard {
        baseline
    } else {
        baseline
            .iter()
            .zip(envelope.iter())
            .map(|(base, env)| base.clamp(env - lower_guard, env + upper_guard))
            .collect::<Vec<_>>()
    };

    Ok(values
        .iter()
        .zip(guarded.iter())
        .map(|(value, base)| (value - base).max(0.0))
        .collect())
}

pub fn baseline_correct_quantile_nonnegative(
    values: &[f64],
    bin_size: usize,
    quantile: f64,
) -> Vec<f64> {
    let baseline = rolling_quantile_baseline(values, bin_size, quantile);
    values
        .iter()
        .zip(baseline.iter())
        .map(|(value, base)| (value - base).max(0.0))
        .collect()
}

fn rolling_quantile_baseline(values: &[f64], bin_size: usize, quantile: f64) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let effective_bin_size = bin_size.max(20);
    let n = values.len();
    let n_bins = (n + effective_bin_size - 1) / effective_bin_size;
    let mut centers = Vec::with_capacity(n_bins);
    let mut base_vals = Vec::with_capacity(n_bins);

    for bin_index in 0..n_bins {
        let start = bin_index * effective_bin_size;
        let end = (start + effective_bin_size).min(n);
        let segment = &values[start..end];
        if segment.is_empty() {
            continue;
        }
        centers.push((start + end - 1) as f64 * 0.5);
        base_vals.push(slice_quantile(segment, quantile));
    }

    if centers.is_empty() {
        return vec![0.0; n];
    }

    let mut baseline = vec![base_vals[0]; n];
    let mut segment_index = 0usize;
    for index in 0..n {
        let x = index as f64;
        while segment_index + 1 < centers.len() && x > centers[segment_index + 1] {
            segment_index += 1;
        }
        if segment_index + 1 >= centers.len() {
            baseline[index] = *base_vals.last().unwrap_or(&base_vals[0]);
            continue;
        }
        let left_x = centers[segment_index];
        let right_x = centers[segment_index + 1];
        let left_y = base_vals[segment_index];
        let right_y = base_vals[segment_index + 1];
        let t = if right_x > left_x {
            (x - left_x) / (right_x - left_x)
        } else {
            0.0
        };
        baseline[index] = left_y + t * (right_y - left_y);
    }

    baseline
}

fn slice_quantile(values: &[f64], quantile: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|left, right| left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal));
    if sorted.len() == 1 {
        return sorted[0];
    }
    let q = quantile.clamp(0.0, 1.0);
    let position = q * (sorted.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        return sorted[lower];
    }
    let weight = position - lower as f64;
    sorted[lower] + weight * (sorted[upper] - sorted[lower])
}

fn stddev(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| {
            let delta = value - mean;
            delta * delta
        })
        .sum::<f64>()
        / values.len() as f64;
    variance.sqrt()
}

fn difference_penalty_bands(
    n: usize,
    lambda: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut lower2 = vec![0.0; n];
    let mut lower1 = vec![0.0; n];
    let mut main = vec![0.0; n];
    let mut upper1 = vec![0.0; n];
    let mut upper2 = vec![0.0; n];

    for index in 0..n {
        main[index] = if index == 0 || index == n - 1 {
            lambda
        } else if index == 1 || index == n - 2 {
            5.0 * lambda
        } else {
            6.0 * lambda
        };
        if index + 1 < n {
            upper1[index] = if index == 0 || index + 1 == n - 1 {
                -2.0 * lambda
            } else {
                -4.0 * lambda
            };
        }
        if index + 2 < n {
            upper2[index] = lambda;
        }
        if index >= 1 {
            lower1[index] = upper1[index - 1];
        }
        if index >= 2 {
            lower2[index] = upper2[index - 2];
        }
    }

    (lower2, lower1, main, upper1, upper2)
}

fn solve_pentadiagonal(
    lower2: &[f64],
    lower1: &[f64],
    main: &[f64],
    upper1: &[f64],
    upper2: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, EngineError> {
    let n = main.len();
    if lower2.len() != n
        || lower1.len() != n
        || upper1.len() != n
        || upper2.len() != n
        || rhs.len() != n
    {
        return Err(EngineError::SignalMath(
            "pentadiagonal system has inconsistent dimensions".to_owned(),
        ));
    }

    let mut alpha = vec![0.0; n];
    let mut gamma = vec![0.0; n];
    let mut delta = vec![0.0; n];
    let mut z = vec![0.0; n];

    alpha[0] = main[0];
    if alpha[0].abs() <= f64::EPSILON {
        return Err(EngineError::SignalMath(
            "singular pentadiagonal system at row 0".to_owned(),
        ));
    }
    if n > 1 {
        gamma[0] = upper1[0] / alpha[0];
    }
    if n > 2 {
        delta[0] = upper2[0] / alpha[0];
    }
    z[0] = rhs[0] / alpha[0];

    if n > 1 {
        alpha[1] = main[1] - lower1[1] * gamma[0];
        if alpha[1].abs() <= f64::EPSILON {
            return Err(EngineError::SignalMath(
                "singular pentadiagonal system at row 1".to_owned(),
            ));
        }
        if n > 2 {
            gamma[1] = (upper1[1] - lower1[1] * delta[0]) / alpha[1];
        }
        if n > 3 {
            delta[1] = upper2[1] / alpha[1];
        }
        z[1] = (rhs[1] - lower1[1] * z[0]) / alpha[1];
    }

    for index in 2..n {
        alpha[index] =
            main[index] - lower2[index] * delta[index - 2] - lower1[index] * gamma[index - 1];
        if alpha[index].abs() <= f64::EPSILON {
            return Err(EngineError::SignalMath(format!(
                "singular pentadiagonal system at row {index}"
            )));
        }
        if index < n - 1 {
            gamma[index] = (upper1[index] - lower1[index] * delta[index - 1]) / alpha[index];
        }
        if index < n - 2 {
            delta[index] = upper2[index] / alpha[index];
        }
        z[index] = (rhs[index] - lower2[index] * z[index - 2] - lower1[index] * z[index - 1])
            / alpha[index];
    }

    let mut solution = vec![0.0; n];
    solution[n - 1] = z[n - 1];
    if n > 1 {
        solution[n - 2] = z[n - 2] - gamma[n - 2] * solution[n - 1];
    }
    if n > 2 {
        for index in (0..=n - 3).rev() {
            solution[index] =
                z[index] - gamma[index] * solution[index + 1] - delta[index] * solution[index + 2];
        }
    }
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::{baseline_correct_nonnegative, find_peaks, peak_score};

    #[test]
    fn peak_finder_respects_height_and_distance() {
        let values = vec![0.0, 1.0, 8.0, 2.0, 0.0, 0.0, 7.5, 1.0, 0.0, 9.0, 0.0];
        let peaks = find_peaks(&values, 5.0, 3);
        let indices = peaks.iter().map(|peak| peak.index).collect::<Vec<_>>();
        assert_eq!(indices, vec![2, 6, 9]);
    }

    #[test]
    fn baseline_correction_preserves_peak_and_clamps_negative_values() {
        let values = (0..120)
            .map(|idx| {
                let drift = 5.0 + idx as f64 * 0.03;
                let peak = if (55..=60).contains(&idx) { 40.0 } else { 0.0 };
                drift + peak
            })
            .collect::<Vec<_>>();
        let corrected = baseline_correct_nonnegative(&values, 0.01, 1_000.0, 50)
            .expect("baseline correction should succeed");
        let max_value = corrected.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min_value = corrected.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(max_value > 10.0);
        assert!(min_value >= 0.0);
    }

    #[test]
    fn peak_finder_prefers_clean_peak_over_blob_shoulder() {
        let shoulder_score = peak_score(520.0, 65.0, 18.0, 430.0);
        let clean_score = peak_score(240.0, 210.0, 4.0, 20.0);
        assert!(clean_score > shoulder_score);
    }
}
