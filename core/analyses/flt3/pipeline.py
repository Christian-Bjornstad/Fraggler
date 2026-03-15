from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from fraggler.fraggler import (
    FsaFile,
    calculate_best_combination_of_size_standard_peaks,
    find_size_standard_peaks,
    fit_size_standard_to_ladder,
    generate_combinations,
    print_green,
    print_warning,
    return_maxium_allowed_distance_between_size_standard_peaks,
)

from core.analysis import (
    MIN_R2_QUALITY,
    _refine_polynomial,
    _select_best_ladder_candidate,
    analyse_fsa_rox,
    compute_ladder_qc_metrics,
    estimate_running_baseline,
)
from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import ASSAY_CONFIG, BP_CORRECTION_OFFSETS, PREFERRED_INJECTION_TIME
from core.analyses.shared_pipeline import finalize_pipeline_run, normalize_pipeline_paths, scan_fsa_files
from core.html_reports import extract_dit_from_name


FLT3_LADDER_QC_THRESHOLD = 0.99
RELEVANT_PEAK_LABELS = {"WT", "MUT", "ITD"}


def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    """Scan recursively for FLT3 .fsa files, excluding water/Vann files."""
    return scan_fsa_files(fsa_dir, mode=mode, recursive=True)


def _preferred_injection_time(meta: dict) -> int:
    assay = meta.get("assay")
    analysis_type = meta.get("analysis_type")
    if analysis_type == "ratio_quant":
        return 1
    if assay == "FLT3-D835":
        return 3
    if analysis_type in PREFERRED_INJECTION_TIME:
        return int(PREFERRED_INJECTION_TIME[analysis_type])
    return int(PREFERRED_INJECTION_TIME.get(assay, meta.get("protocol_injection_time", meta.get("injection_time", 0))))


def _candidate_sort_key(item: tuple[Path, dict], preferred_injection: int) -> tuple[int, int, str, str]:
    path, meta = item
    injection = int(meta.get("injection_time", 0) or 0)
    return (
        0 if injection == preferred_injection else 1,
        abs(injection - preferred_injection),
        meta.get("source_run_dir", ""),
        path.name,
    )


def _calculate_auc(trace: np.ndarray, time_idx: np.ndarray) -> float:
    if time_idx.size == 0:
        return 0.0
    clipped = time_idx[(time_idx >= 0) & (time_idx < trace.size)]
    if clipped.size == 0:
        return 0.0
    return float(trace[clipped].sum())


def _peak_area_half_width_bp(assay: str, label: str, center_bp: float) -> float:
    if assay == "FLT3-D835":
        if label == "WT":
            return 1.2
        if label == "MUT":
            return 0.5
        if abs(center_bp - 150.0) <= 6.0:
            return 0.8
        return 0.8
    return 5.0


def _calculate_peak_area(
    trace: np.ndarray,
    time_all: np.ndarray,
    bp_all: np.ndarray,
    center_bp: float,
    assay: str,
    label: str,
) -> float:
    half_width_bp = _peak_area_half_width_bp(assay, label, center_bp)
    auc_mask = (bp_all >= center_bp - half_width_bp) & (bp_all <= center_bp + half_width_bp)
    return _calculate_auc(trace, time_all[auc_mask])


def _assay_positive_ratio(assay: str) -> float:
    return float(ASSAY_CONFIG.get(assay, {}).get("positive_ratio", 0.01))


def _bp_in_ranges(bp: float, ranges: list[tuple[float, float]] | None) -> bool:
    if not ranges:
        return False
    return any(start <= bp <= end for start, end in ranges)


def _apply_bp_offset(fsa: FsaFile, assay: str) -> None:
    offset = float(BP_CORRECTION_OFFSETS.get(assay, 0.0))
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if not offset or sample_data is None or sample_data.empty:
        return
    adjusted = sample_data.copy()
    adjusted["basepairs"] = adjusted["basepairs"].astype(float) + offset
    fsa.sample_data_with_basepairs = adjusted


def _infer_sizing_method(fsa: FsaFile) -> str:
    if hasattr(fsa, "_flt3_sizing_method"):
        return str(getattr(fsa, "_flt3_sizing_method"))
    model = getattr(fsa, "ladder_model", None)
    model_name = type(model).__name__
    if model_name == "Pipeline":
        return "spline"
    if model_name == "LinearRegression":
        return "polynomial_refinement"
    return "unknown"


def _detect_peaks(
    fsa: FsaFile,
    assay: str,
    wt_bp: float,
    trace: np.ndarray,
    mut_bp: float | None = None,
    analysis_type: str | None = None,
) -> pd.DataFrame:
    """Detect WT and mutant peaks and estimate their corrected AUC."""
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if sample_data is None or sample_data.empty:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    time_all = sample_data["time"].astype(int).to_numpy()
    bp_all = sample_data["basepairs"].to_numpy()

    bp_min, bp_max = 50.0, 1000.0
    peaks: list[dict] = []
    assay_cfg = ASSAY_CONFIG.get(assay, {})

    from scipy.signal import find_peaks

    mask = (bp_all >= bp_min) & (bp_all <= bp_max)
    if not mask.any():
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    valid_time = time_all[mask]
    valid_time = valid_time[(valid_time >= 0) & (valid_time < trace.size)]
    if valid_time.size < 3:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    y_win = trace[time_all[mask]]
    bp_win = bp_all[mask]

    peak_height_min = float(assay_cfg.get("peak_height_min", 200))
    peak_distance = int(assay_cfg.get("peak_distance", 20))
    peak_idx, _ = find_peaks(y_win, height=peak_height_min, distance=peak_distance)

    wt_tol = 4.0 if assay == "FLT3-ITD" else 2.0 if assay in {"FLT3-D835", "NPM1"} else 5.0
    itd_min_bp = float(ASSAY_CONFIG.get("FLT3-ITD", {}).get("itd_min_bp", wt_bp + 4.9)) if wt_bp else None
    wt_range = assay_cfg.get("wt_range")
    mut_ranges = assay_cfg.get("mut_ranges")

    channel_traces: dict[str, np.ndarray] = {}
    for ch in ["DATA1", "DATA2", "DATA3"]:
        if ch not in fsa.fsa:
            continue
        raw = np.asarray(fsa.fsa[ch]).astype(float)
        baseline = estimate_running_baseline(raw, bin_size=5000, quantile=0.01)
        channel_traces[ch] = np.maximum(raw - baseline, 0.0)

    for idx in peak_idx:
        p_bp = float(bp_win[idx])
        p_h = float(y_win[idx])

        label = "unspecific"
        if _bp_in_ranges(p_bp, [wt_range] if wt_range else None) or (wt_bp and abs(p_bp - wt_bp) < wt_tol):
            label = "WT"
        elif assay == "FLT3-ITD" and itd_min_bp is not None and p_bp >= itd_min_bp:
            label = "ITD"
        elif _bp_in_ranges(p_bp, mut_ranges) or (mut_bp and abs(p_bp - mut_bp) < (wt_tol + 2)):
            label = "MUT"

        p_area = _calculate_peak_area(
            trace=trace,
            time_all=time_all,
            bp_all=bp_all,
            center_bp=p_bp,
            assay=assay,
            label=label,
        )

        peak_info = {
            "basepairs": p_bp,
            "peaks": p_h,
            "area": p_area,
            "label": label,
            "keep": True,
        }
        for ch, corr_trace in channel_traces.items():
            peak_info[f"area_{ch}"] = _calculate_peak_area(
                trace=corr_trace,
                time_all=time_all,
                bp_all=bp_all,
                center_bp=p_bp,
                assay=assay,
                label=label,
            )
        peaks.append(peak_info)

    if not peaks:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
    return pd.DataFrame(peaks)


def _combine_peak_traces(fsa: FsaFile, peak_channels: list[str], primary_channel: str) -> np.ndarray:
    combined_trace = None
    for ch in peak_channels:
        if ch not in fsa.fsa:
            continue
        raw = np.asarray(fsa.fsa[ch]).astype(float)
        baseline = estimate_running_baseline(raw, bin_size=5000, quantile=0.01)
        corrected = np.maximum(raw - baseline, 0.0)
        combined_trace = corrected if combined_trace is None else combined_trace + corrected

    if combined_trace is not None:
        return combined_trace

    if primary_channel in fsa.fsa:
        return np.asarray(fsa.fsa[primary_channel]).astype(float)

    first_trace = next(iter(fsa.fsa.values()))
    return np.zeros(len(first_trace), dtype=float)


def _peak_qc_status(peaks: pd.DataFrame, group: str) -> tuple[bool, str]:
    if group == "negative_control":
        return True, "negative_control"
    if peaks.empty:
        return False, "no_peaks"
    relevant = peaks[peaks.label.isin(RELEVANT_PEAK_LABELS)]
    if relevant.empty:
        return False, "no_relevant_peaks"
    return True, "ok"


def _attempt_lenient_rox_fit(fsa_path: Path, sample_channel: str) -> FsaFile | None:
    configs = [{"min_h": 5, "min_d": 3}]
    for cfg in configs:
        try:
            fsa = FsaFile(
                file=str(fsa_path),
                ladder="GS500ROX",
                sample_channel=sample_channel,
                min_distance_between_peaks=cfg["min_d"],
                min_size_standard_height=cfg["min_h"],
                size_standard_channel="DATA4",
            )
            fsa = find_size_standard_peaks(fsa)
            ss_peaks = getattr(fsa, "size_standard_peaks", None)
            if ss_peaks is None or getattr(ss_peaks, "shape", [0])[0] < 3:
                continue

            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=1.5)
            for _ in range(20):
                fsa = generate_combinations(fsa)
                best = getattr(fsa, "best_size_standard_combinations", None)
                if best is not None and best.shape[0] > 0:
                    break
                fsa.maxium_allowed_distance_between_size_standard_peaks += 10

            best = getattr(fsa, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                continue

            selected_fit = _select_best_ladder_candidate(fsa)
            if selected_fit is not None:
                fsa = selected_fit
            else:
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)

            if not getattr(fsa, "fitted_to_model", False):
                continue

            qc = compute_ladder_qc_metrics(fsa)
            if qc["r2"] >= MIN_R2_QUALITY:
                setattr(fsa, "_flt3_sizing_method", "spline_lenient")
                return fsa

            refined = _refine_polynomial(fsa, "ROX-lenient", fsa_path)
            if refined is not None:
                setattr(refined, "_flt3_sizing_method", "polynomial_refinement_lenient")
                return refined

            setattr(fsa, "_flt3_sizing_method", "spline_lenient")
            return fsa
        except Exception:
            continue

    return None


def _analyse_fsa_candidate(fsa_path: Path, sample_channel: str, assay: str) -> FsaFile | None:
    fsa = analyse_fsa_rox(fsa_path, sample_channel)
    if fsa is not None:
        setattr(fsa, "_flt3_sizing_method", _infer_sizing_method(fsa))
        return fsa
    if assay == "FLT3-D835":
        return _attempt_lenient_rox_fit(fsa_path, sample_channel)
    return None


def _build_entry_from_candidate(fsa_path: Path, meta: dict) -> dict | None:
    fsa = _analyse_fsa_candidate(fsa_path, meta["primary_peak_channel"], meta["assay"])
    if fsa is None:
        return None

    _apply_bp_offset(fsa, meta["assay"])
    combined_trace = _combine_peak_traces(
        fsa=fsa,
        peak_channels=meta.get("peak_channels", [meta["primary_peak_channel"]]),
        primary_channel=meta["primary_peak_channel"],
    )
    peaks = _detect_peaks(
        fsa=fsa,
        assay=meta["assay"],
        wt_bp=meta["wt_bp"],
        trace=combined_trace,
        mut_bp=meta.get("mut_bp"),
        analysis_type=meta.get("analysis_type"),
    )

    metrics = compute_ladder_qc_metrics(fsa)
    ladder_qc_status = "ok" if float(metrics.get("r2", float("nan"))) > FLT3_LADDER_QC_THRESHOLD else "ladder_qc_failed"
    peak_qc_pass, peak_qc_reason = _peak_qc_status(peaks, meta.get("group", "sample"))

    return {
        "fsa": fsa,
        "peaks_by_channel": {meta["primary_peak_channel"]: peaks},
        "trace_channels": meta["trace_channels"],
        "primary_peak_channel": meta["primary_peak_channel"],
        "ymax": float(np.max(combined_trace)) * 1.1 if combined_trace.size and np.any(combined_trace) else 1000.0,
        "assay": meta["assay"],
        "analysis_type": meta["analysis_type"],
        "parallel": meta.get("parallel"),
        "well_id": meta.get("well_id"),
        "specimen_id": meta.get("specimen_id"),
        "selection_key": meta.get("selection_key"),
        "group": meta["group"],
        "ladder": fsa.ladder,
        "bp_min": meta["bp_min"],
        "bp_max": meta["bp_max"],
        "dit": extract_dit_from_name(fsa.file_name),
        "ladder_qc_status": ladder_qc_status,
        "ladder_r2": float(metrics.get("r2", np.nan)),
        "n_ladder_steps": metrics.get("n_ladder_steps"),
        "n_size_standard_peaks": metrics.get("n_size_standard_peaks"),
        "injection_time": meta["injection_time"],
        "selected_injection": f"{int(meta['injection_time'])}s",
        "selected_injection_time": int(meta["injection_time"]),
        "preferred_injection_time": _preferred_injection_time(meta),
        "protocol_injection_time": meta.get("protocol_injection_time", meta["injection_time"]),
        "source_run_dir": meta.get("source_run_dir", ""),
        "run_name": meta.get("run_name", ""),
        "run_date": meta.get("run_date", ""),
        "run_time": meta.get("run_time", ""),
        "injection_protocol": meta.get("injection_protocol", ""),
        "selection_reason": "",
        "alternate_injections": [],
        "alternate_injections_summary": "",
        "sizing_method": _infer_sizing_method(fsa),
        "peak_qc_pass": peak_qc_pass,
        "peak_qc_status": peak_qc_reason,
    }


def _candidate_audit_record(path: Path, meta: dict, status: str, reason: str) -> dict:
    return {
        "file": path.name,
        "injection_time": int(meta.get("injection_time", 0) or 0),
        "selected_injection": f"{int(meta.get('injection_time', 0) or 0)}s",
        "source_run_dir": meta.get("source_run_dir", ""),
        "status": status,
        "reason": reason,
    }


def _select_best_entry(candidates: list[tuple[Path, dict]]) -> dict | None:
    if not candidates:
        return None

    preferred_injection = _preferred_injection_time(candidates[0][1])
    ordered = sorted(candidates, key=lambda item: _candidate_sort_key(item, preferred_injection))

    audit_records: list[dict] = []
    best_available: tuple[dict, str] | None = None
    preferred_reason = "preferred injection unavailable"

    for index, (path, meta) in enumerate(ordered):
        same_as_preferred = int(meta.get("injection_time", 0) or 0) == preferred_injection
        entry = _build_entry_from_candidate(path, meta)
        if entry is None:
            reason = "ladder_fit_failed"
            if same_as_preferred:
                preferred_reason = reason
            audit_records.append(_candidate_audit_record(path, meta, "rejected", reason))
            continue

        acceptable = entry["ladder_qc_status"] == "ok" and entry["peak_qc_pass"]
        candidate_reason = "qc_pass"
        if not acceptable:
            if entry["ladder_qc_status"] != "ok":
                candidate_reason = entry["ladder_qc_status"]
            else:
                candidate_reason = entry["peak_qc_status"]
            if same_as_preferred:
                preferred_reason = candidate_reason
            if best_available is None:
                best_available = (entry, candidate_reason)
            audit_records.append(_candidate_audit_record(path, meta, "rejected", candidate_reason))
            continue

        if same_as_preferred:
            selection_reason = f"Preferred {preferred_injection}s injection selected"
        elif any(int(m.get("injection_time", 0) or 0) == preferred_injection for _, m in candidates):
            selection_reason = f"Preferred {preferred_injection}s failed ({preferred_reason}); selected {entry['selected_injection']} fallback"
        else:
            selection_reason = f"Preferred {preferred_injection}s unavailable; selected {entry['selected_injection']}"

        entry["selection_reason"] = selection_reason
        entry["preferred_injection_time"] = preferred_injection
        entry["alternate_injections"] = [
            record for record in audit_records
            if record["file"] != path.name
        ] + [
            _candidate_audit_record(other_path, other_meta, "not_selected", "preferred candidate passed")
            for other_path, other_meta in ordered[index + 1:]
        ]
        entry["alternate_injections_summary"] = "; ".join(
            f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
            for alt in entry["alternate_injections"]
        )
        return entry

    if best_available is None:
        return None

    entry, best_reason = best_available
    entry["selection_reason"] = f"No candidate passed QC; kept {entry['selected_injection']} ({best_reason})"
    entry["preferred_injection_time"] = preferred_injection
    entry["alternate_injections"] = audit_records
    entry["alternate_injections_summary"] = "; ".join(
        f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
        for alt in audit_records
        if alt["file"] != entry["fsa"].file_name
    )
    return entry


def _calculate_ratios(entries: list[dict]) -> None:
    """Calculate FLT3 mutant ratios and store explicit numerator/denominator fields."""
    for entry in entries:
        peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
        if peaks.empty:
            entry["ratio"] = 0.0
            entry["ratio_numerator_area"] = 0.0
            entry["ratio_denominator_area"] = 0.0
            entry["mutant_fraction"] = 0.0
            continue

        wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False)
        mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])]

        if wt_rows.empty:
            entry["ratio"] = 0.0
            entry["ratio_numerator_area"] = 0.0
            entry["ratio_denominator_area"] = 0.0
            entry["mutant_fraction"] = 0.0
            continue

        wt_main = wt_rows.iloc[0]
        wt_area = float(wt_main.area)

        if entry["assay"] == "FLT3-ITD":
            mut_area = float(mut_rows.area.sum()) if not mut_rows.empty else 0.0
        elif not mut_rows.empty:
            mut_area = float(mut_rows.sort_values("area", ascending=False).iloc[0].area)
        else:
            mut_area = 0.0

        entry["ratio_numerator_area"] = mut_area
        entry["ratio_denominator_area"] = wt_area
        entry["ratio"] = (mut_area / wt_area) if wt_area > 0 else 0.0
        entry["mutant_fraction"] = (mut_area / (mut_area + wt_area)) if (mut_area + wt_area) > 0 else 0.0


def _summarize_peak_areas(entry: dict) -> tuple[float, float]:
    return float(entry.get("ratio_denominator_area", 0.0)), float(entry.get("ratio_numerator_area", 0.0))


def _summarize_detected_peaks(entry: dict) -> dict:
    peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
    wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False) if not peaks.empty else pd.DataFrame()
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].sort_values("basepairs") if not peaks.empty else pd.DataFrame()

    wt_bp = float(wt_rows.iloc[0].basepairs) if not wt_rows.empty else np.nan
    wt_area = float(wt_rows.iloc[0].area) if not wt_rows.empty else 0.0
    mut_bps = [round(float(v), 2) for v in mut_rows.basepairs.tolist()] if not mut_rows.empty else []
    mut_areas = [round(float(v), 2) for v in mut_rows.area.tolist()] if not mut_rows.empty else []
    mut_main_bp = float(mut_rows.sort_values("area", ascending=False).iloc[0].basepairs) if not mut_rows.empty else np.nan
    mut_main_area = float(mut_rows.sort_values("area", ascending=False).iloc[0].area) if not mut_rows.empty else 0.0

    return {
        "wt_bp": wt_bp,
        "wt_area": wt_area,
        "mut_bps": mut_bps,
        "mut_areas": mut_areas,
        "mut_area_total": float(sum(mut_areas)),
        "mut_main_bp": mut_main_bp,
        "mut_main_area": mut_main_area,
    }


def _interpret_entry(entry: dict) -> str:
    assay = entry["assay"]
    ratio = float(entry.get("ratio", 0.0))
    peak_summary = _summarize_detected_peaks(entry)
    positive_ratio = _assay_positive_ratio(assay)

    if assay == "FLT3-ITD":
        if ratio >= positive_ratio:
            return "Positiv FLT3-ITD"
        if peak_summary["mut_bps"]:
            return "Negativ FLT3-ITD - lavniva dokumentert"
        return "Ingen FLT3-ITD pavist"
    if assay == "FLT3-D835":
        if ratio >= positive_ratio:
            return "Positiv FLT3-D835"
        if peak_summary["mut_bps"]:
            return "FLT3-D835 under positiv grense - dokumentert"
        return "Ingen FLT3-D835 pavist"
    if assay == "NPM1":
        if ratio >= positive_ratio:
            return "Positiv NPM1"
        if peak_summary["mut_bps"]:
            return "Mulig NPM1 - vurder manuelt"
        return "Ingen NPM1-mutasjon pavist"
    return "Ingen tolkning"


def generate_flt3_peak_report(entries: list[dict], outdir: Path) -> None:
    rows = []
    for entry in entries:
        peak_summary = _summarize_detected_peaks(entry)
        rows.append(
            {
                "DIT": entry.get("dit") or "",
                "File": entry["fsa"].file_name,
                "Assay": entry["assay"],
                "Group": entry.get("group") or "",
                "Parallel": entry.get("parallel") or "",
                "Well": entry.get("well_id") or "",
                "Treatment": entry.get("analysis_type") or "",
                "SelectedInjection": entry.get("selected_injection") or "",
                "PreferredInjection": f"{int(entry.get('preferred_injection_time', 0) or 0)}s" if entry.get("preferred_injection_time") else "",
                "SelectionReason": entry.get("selection_reason") or "",
                "SourceRunDir": entry.get("source_run_dir") or "",
                "AlternateInjections": entry.get("alternate_injections_summary") or "",
                "SizingMethod": entry.get("sizing_method") or "",
                "InjectionTime": entry.get("injection_time"),
                "WT_bp": round(float(peak_summary["wt_bp"]), 2) if not np.isnan(peak_summary["wt_bp"]) else "",
                "WT_Area": round(peak_summary["wt_area"], 2),
                "Mutant_bp": ", ".join(f"{bp:.2f}" for bp in peak_summary["mut_bps"]),
                "Mutant_Area": ", ".join(f"{area:.2f}" for area in peak_summary["mut_areas"]),
                "Mutant_Area_Total": round(peak_summary["mut_area_total"], 2),
                "RatioNumeratorArea": round(float(entry.get("ratio_numerator_area", 0.0)), 2),
                "RatioDenominatorArea": round(float(entry.get("ratio_denominator_area", 0.0)), 2),
                "Ratio": round(float(entry.get("ratio", 0.0)), 4),
                "MutantFractionMutPlusWT": round(float(entry.get("mutant_fraction", 0.0)), 4),
                "Interpretation": _interpret_entry(entry),
                "LadderQC": entry.get("ladder_qc_status", ""),
                "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
            }
        )

    if not rows:
        return

    df = pd.DataFrame(rows)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "Final_Detailed_Peak_Report.csv"
    df.to_csv(csv_path, index=False)
    print_green(f"FLT3 detailed peak report saved to {csv_path}")


def generate_flt3_qc_report(entries: list[dict], outdir: Path) -> None:
    qc_data = []

    for entry in entries:
        group = entry.get("group")
        if group not in ["negative_control", "positive_control", "reactive_control"]:
            continue

        peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
        wt_peaks = peaks[peaks.label == "WT"]
        mut_peaks = peaks[peaks.label.isin(["MUT", "ITD"])]
        wt_area, mut_area = _summarize_peak_areas(entry)
        ratio = float(entry.get("ratio", 0.0))
        assay_cfg = ASSAY_CONFIG.get(entry["assay"], {})
        min_wt_area = float(assay_cfg.get("control_wt_min_area", 0.0))
        min_ratio = _assay_positive_ratio(entry["assay"])

        status = "FAIL"
        details = ""

        if group == "negative_control":
            if mut_peaks.empty:
                status = "PASS"
            else:
                details = f"Unexpected mutant peaks found: {mut_peaks.basepairs.tolist()}"
        elif group == "reactive_control":
            if not wt_peaks.empty and wt_area >= min_wt_area:
                status = "PASS"
            elif wt_peaks.empty:
                details = "No WT peak detected"
            else:
                details = f"WT area below threshold ({wt_area:.0f} < {min_wt_area:.0f})"
        elif group == "positive_control":
            if not mut_peaks.empty and (ratio >= min_ratio or wt_peaks.empty):
                status = "PASS"
            elif mut_peaks.empty:
                details = "No mutant/ITD peak detected"
            else:
                details = f"Mutant ratio below threshold ({ratio:.4f} < {min_ratio:.4f})"

        qc_data.append(
            {
                "File": entry["fsa"].file_name,
                "ControlGroup": group,
                "Assay": entry["assay"],
                "Well": entry.get("well_id") or "",
                "SelectedInjection": entry.get("selected_injection") or "",
                "SelectionReason": entry.get("selection_reason") or "",
                "InjectionTime": entry.get("injection_time"),
                "WT_Area": round(wt_area, 2),
                "Mutant_Area": round(mut_area, 2),
                "Ratio": round(ratio, 4),
                "Status": status,
                "Details": details,
            }
        )

    if not qc_data:
        return

    df = pd.DataFrame(qc_data)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "QC_FLT3_Injections.csv"
    df.to_csv(csv_path, index=False)
    print_green(f"FLT3 QC report saved to {csv_path}")


def generate_flt3_bp_validation_report(entries: list[dict], outdir: Path) -> None:
    rows = []
    for entry in entries:
        peak_summary = _summarize_detected_peaks(entry)
        assay_cfg = ASSAY_CONFIG.get(entry["assay"], {})
        if not np.isnan(peak_summary["wt_bp"]):
            expected_wt = float(assay_cfg.get("wt_bp", np.nan))
            rows.append(
                {
                    "DIT": entry.get("dit") or "",
                    "File": entry["fsa"].file_name,
                    "Assay": entry["assay"],
                    "Group": entry.get("group") or "",
                    "Well": entry.get("well_id") or "",
                    "InjectionTime": entry.get("injection_time"),
                    "SelectedInjection": entry.get("selected_injection") or "",
                    "SizingMethod": entry.get("sizing_method") or "",
                    "PeakType": "WT",
                    "ExpectedBP": round(expected_wt, 2) if np.isfinite(expected_wt) else "",
                    "ObservedBP": round(float(peak_summary["wt_bp"]), 2),
                    "DeltaBP": round(float(peak_summary["wt_bp"]) - expected_wt, 2) if np.isfinite(expected_wt) else "",
                    "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
                }
            )

        if entry["assay"] in {"FLT3-D835", "NPM1"} and not np.isnan(peak_summary["mut_main_bp"]):
            expected_mut = float(assay_cfg.get("mut_bp", np.nan))
            rows.append(
                {
                    "DIT": entry.get("dit") or "",
                    "File": entry["fsa"].file_name,
                    "Assay": entry["assay"],
                    "Group": entry.get("group") or "",
                    "Well": entry.get("well_id") or "",
                    "InjectionTime": entry.get("injection_time"),
                    "SelectedInjection": entry.get("selected_injection") or "",
                    "SizingMethod": entry.get("sizing_method") or "",
                    "PeakType": "MUT",
                    "ExpectedBP": round(expected_mut, 2) if np.isfinite(expected_mut) else "",
                    "ObservedBP": round(float(peak_summary["mut_main_bp"]), 2),
                    "DeltaBP": round(float(peak_summary["mut_main_bp"]) - expected_mut, 2) if np.isfinite(expected_mut) else "",
                    "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
                }
            )

    if not rows:
        return

    df = pd.DataFrame(rows)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "FLT3_BP_Validation.csv"
    df.to_csv(csv_path, index=False)
    print_green(f"FLT3 bp validation report saved to {csv_path}")


def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> list[dict] | None:
    """
    Kjor FLT3-pipeline pa alle .fsa-filer i fsa_dir.
    """
    fsa_dir, assay_dir = normalize_pipeline_paths(fsa_dir, base_outdir, assay_folder_name)

    raw_files = _scan_files(fsa_dir, mode=mode)
    classified = []
    for path in raw_files:
        meta = classify_fsa(path)
        if meta:
            classified.append((path, meta))

    if not classified:
        return [] if return_entries else None

    groups: dict[str, list[tuple[Path, dict]]] = defaultdict(list)
    for path, meta in classified:
        groups[meta["selection_key"]].append((path, meta))

    entries = []
    for selection_key, candidates in sorted(groups.items()):
        entry = _select_best_entry(candidates)
        if entry is None:
            first_file = candidates[0][0].name if candidates else selection_key
            print_warning(f"FLT3 selection failed for {first_file}")
            continue
        entries.append(entry)

    if not entries:
        return [] if return_entries else None

    _calculate_ratios(entries)
    generate_flt3_qc_report(entries, assay_dir)
    generate_flt3_peak_report(entries, assay_dir)
    generate_flt3_bp_validation_report(entries, assay_dir)

    return finalize_pipeline_run(
        entries,
        assay_dir,
        return_entries=return_entries,
        make_dit_reports=make_dit_reports,
        mode=mode,
    )
