use serde::{Deserialize, Serialize};

use crate::engine::EngineError;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Peak {
    pub index: usize,
    pub height: f64,
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
            candidates.push(Peak {
                index,
                height: current,
            });
        }
    }

    candidates.sort_by(|left, right| {
        right
            .height
            .partial_cmp(&left.height)
            .unwrap_or(std::cmp::Ordering::Equal)
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
    use super::{baseline_correct_nonnegative, find_peaks};

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
}
