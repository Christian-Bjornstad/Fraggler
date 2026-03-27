"""
Fraggler QC — Marker tracking, filename parsing, peak finding utilities.
"""
from __future__ import annotations

import re
from pathlib import Path
from html import escape as html_escape_builtin

import numpy as np

from core.qc.qc_rules import QCRules, normalize_assay_qc
from core.analysis import estimate_running_baseline
from core.utils import strip_stage_prefix, CONTROL_PREFIX_RE


def parse_run_code_from_filename(name: str) -> str | None:
    """
    Henter run_code som siste token før .fsa.
    Eksempel: PK1_TCRgA_120126_E05_H9C0U3SI.fsa -> H9C0U3SI
    """
    if not name:
        return None
    stem = Path(name).stem  # uten .fsa
    parts = stem.split("_")
    if len(parts) < 2:
        return None
    run_code = parts[-1].strip()
    return run_code.upper() if run_code else None


def make_run_key(filename: str) -> str:
    """
    Stabil run-identifikator: dato + run_code.
    Eksempel: 2026-01-12_H9C0U3SI
    """
    d = parse_pcr_date_from_filename(filename)
    c = parse_run_code_from_filename(filename)

    if d and c:
        return f"{d}_{c}"
    if d:
        return d
    if c:
        return c
    return "UNKNOWN"

DATE8_RE = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{4})(?!\d)")  # ddmmyyyy
WELL_RE = re.compile(r"_([A-H]\d{2})_", re.IGNORECASE)        # _G09_
BATCH_RE = re.compile(r"_([A-Z]\d{6}[A-Z])(?:\.fsa)?$", re.IGNORECASE)  # _C991475U.fsa

DATE6_RE = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{2})(?!\d)")      # ddmmyy

def parse_pcr_date_from_filename(name: str) -> str | None:
    """
    Henter dato fra filnavn og returnerer ISO-format YYYY-MM-DD.
    Støtter:
      - ddmmyy  (f.eks. 120126 -> 2026-01-12)
      - ddmmyyyy (f.eks. 06022026 -> 2026-02-06)
    """
    s = name or ""
    s = strip_stage_prefix(s)

    # 1) Prøv ddmmyyyy først (8 siffer)
    m8 = DATE8_RE.search(s)
    if m8:
        dd, mm, yyyy = m8.group(1), m8.group(2), m8.group(3)
        return f"{yyyy}-{mm}-{dd}"

    # 2) Prøv ddmmyy (6 siffer)
    m6 = DATE6_RE.search(s)
    if m6:
        dd, mm, yy = m6.group(1), m6.group(2), m6.group(3)
        # Tolker 00-79 => 2000-2079 (juster hvis dere trenger)
        yyyy = 2000 + int(yy)
        return f"{yyyy:04d}-{mm}-{dd}"

    return None


def parse_well_from_filename(name: str) -> str | None:
    m = WELL_RE.search(name or "")
    return m.group(1).upper() if m else None

def parse_batch_from_filename(name: str) -> str | None:
    # typisk _C991475U.fsa
    m = BATCH_RE.search(name or "")
    return m.group(1).upper() if m else None

def html_escape(s: str) -> str:
    return (s or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")



def control_id_from_filename(filename: str) -> str:
    clean_name = strip_stage_prefix(filename or "")
    m = CONTROL_PREFIX_RE.search(clean_name)
    return m.group(1).upper() if m else "UNKNOWN"


def worst_grade(a: str, b: str) -> str:
    order = {"OK": 0, "WARN": 1, "FAIL": 2, "NA": 1}
    return a if order.get(a, 1) >= order.get(b, 1) else b


def ladder_qc_grade(r2: float | None, rules: QCRules) -> tuple[str, str]:
    if r2 is None or not np.isfinite(r2):
        return ("FAIL", "Ugyldig/manglende R²")

    r2_grade = "OK"
    if r2 >= rules.min_r2_ok:
        r2_grade = "OK"
    elif r2 >= rules.min_r2_warn:
        r2_grade = "WARN"
    else:
        r2_grade = "FAIL"

    return (r2_grade, f"R²={r2:.4f}({r2_grade})")


LOCAL_SCORE_DISTANCE_WEIGHT = 4.0


NON_SL_SAMPLE_TUNING: dict[str, dict[str, float | str]] = {
    "FR1": {"sample_window_bp": 2.0, "fallback_window_bp": 5.0, "selection_mode": "legacy"},
    "FR2": {"sample_window_bp": 2.0, "fallback_window_bp": 5.0, "selection_mode": "legacy"},
    "FR3": {"sample_window_bp": 3.0, "fallback_window_bp": 8.0, "selection_mode": "score"},
    "DHJH_D": {"sample_window_bp": 3.0, "fallback_window_bp": 9.0, "selection_mode": "score"},
    "DHJH_E": {"sample_window_bp": 2.0, "fallback_window_bp": 6.0, "selection_mode": "score"},
    "IGK_279": {"sample_window_bp": 2.0, "fallback_window_bp": 8.0, "selection_mode": "score"},
    "IGK_150": {"sample_window_bp": 2.0, "fallback_window_bp": 8.0, "selection_mode": "score"},
    "KDE": {"sample_window_bp": 4.0, "fallback_window_bp": 12.0, "selection_mode": "score"},
    "TCRbA": {"sample_window_bp": 2.0, "fallback_window_bp": 7.0, "selection_mode": "legacy"},
    "TCRbB": {"sample_window_bp": 2.0, "fallback_window_bp": 7.0, "selection_mode": "legacy"},
    "TCRbC": {"sample_window_bp": 2.0, "fallback_window_bp": 7.0, "selection_mode": "legacy"},
    "TCRgA": {"sample_window_bp": 4.0, "fallback_window_bp": 9.0, "selection_mode": "score"},
    "TCRgB": {"sample_window_bp": 4.0, "fallback_window_bp": 9.0, "selection_mode": "score"},
}

NON_SL_TARGET_TUNING: dict[float, dict[str, float | str]] = {
    325.0: NON_SL_SAMPLE_TUNING["FR1"],
    260.0: NON_SL_SAMPLE_TUNING["FR2"],
    145.0: NON_SL_SAMPLE_TUNING["FR3"],
    139.0: NON_SL_SAMPLE_TUNING["DHJH_D"],
    109.0: NON_SL_SAMPLE_TUNING["DHJH_E"],
    279.0: NON_SL_SAMPLE_TUNING["IGK_279"],
    150.0: NON_SL_SAMPLE_TUNING["IGK_150"],
    287.0: NON_SL_SAMPLE_TUNING["KDE"],
    377.0: NON_SL_SAMPLE_TUNING["KDE"],
    265.0: NON_SL_SAMPLE_TUNING["TCRbA"],
    254.0: NON_SL_SAMPLE_TUNING["TCRbB"],
    311.0: NON_SL_SAMPLE_TUNING["TCRbC"],
    249.0: NON_SL_SAMPLE_TUNING["TCRgA"],
    212.0: NON_SL_SAMPLE_TUNING["TCRgA"],
    163.0: NON_SL_SAMPLE_TUNING["TCRgA"],
    115.0: NON_SL_SAMPLE_TUNING["TCRgB"],
    178.0: NON_SL_SAMPLE_TUNING["TCRgB"],
}


def _merge_marker_tuning(base_window_bp: float, tuning: dict[str, float | str] | None) -> dict[str, float | str]:
    if not tuning:
        return {
            "sample_window_bp": float(base_window_bp),
            "fallback_window_bp": float(max(base_window_bp, 8.0)),
            "selection_mode": "legacy",
        }

    sample_window_bp = float(tuning.get("sample_window_bp", base_window_bp))
    respect_exact_window = str(tuning.get("respect_exact_sample_window", "")).strip().lower() in {"1", "true", "yes", "y"}
    if not respect_exact_window:
        sample_window_bp = max(float(base_window_bp), sample_window_bp)
    fallback_window_bp = max(sample_window_bp, float(tuning.get("fallback_window_bp", sample_window_bp)))
    return {
        "sample_window_bp": sample_window_bp,
        "fallback_window_bp": fallback_window_bp,
        "selection_mode": str(tuning.get("selection_mode", "legacy")),
    }


def _tuning_for_marker_target(target_bp: float, window_bp: float) -> dict[str, float | str]:
    target = round(float(target_bp), 1)
    if target in {100.0, 200.0, 300.0, 400.0, 600.0} and float(window_bp) >= 10.0:
        return {
            "sample_window_bp": float(window_bp),
            "fallback_window_bp": float(max(window_bp, 8.0)),
            "selection_mode": "legacy",
        }
    return _merge_marker_tuning(float(window_bp), NON_SL_TARGET_TUNING.get(target))


def _candidate_score(
    candidate: dict,
    target_bp: float,
    reference_window_bp: float,
    reference_height: float,
) -> float:
    if not candidate.get("ok", False):
        return float("-inf")

    height = float(candidate.get("height", 0.0) or 0.0)
    found_bp = float(candidate.get("found_bp", target_bp) or target_bp)
    delta = abs(found_bp - float(target_bp))
    ref_window = max(float(reference_window_bp), 1.0)
    ref_height = max(float(reference_height), 1.0)
    edge_penalty = 0.35 if delta >= max(ref_window * 0.8, ref_window - 0.5) else 0.0
    return (height / ref_height) * 2.0 - (delta / ref_window) - edge_penalty


def _local_candidate_score(
    candidate: dict,
    target_bp: float,
    reference_window_bp: float,
    reference_height: float,
) -> float:
    if not candidate.get("ok", False):
        return float("-inf")

    height = float(candidate.get("height", 0.0) or 0.0)
    found_bp = float(candidate.get("found_bp", target_bp) or target_bp)
    delta = abs(found_bp - float(target_bp))
    ref_window = max(float(reference_window_bp), 1.0)
    ref_height = max(float(reference_height), 1.0)
    edge_penalty = 0.25 if delta >= max(ref_window * 0.8, ref_window - 0.5) else 0.0
    return (height / ref_height) - (LOCAL_SCORE_DISTANCE_WEIGHT * (delta / ref_window)) - edge_penalty


# ======================================================================
# MARKER KONFIG (forventede peaks + ladder bp som skal trackes)
# - "kind": "sample" eller "ladder"
# - "expected_bp": nominal bp vi leter rundt
# - "channel": "primary" eller f.eks. "DATA4"/"DATA105"
# - "window_bp": bp-halvvindu for søk
# ======================================================================

def markers_for_entry(entry: dict, rules: QCRules) -> list[dict]:
    """
    Returnerer markører som skal trackes for QC.
    Kun PK/PK1/PK2 får markører.
    """
    fsa = entry["fsa"]
    ctrl = control_id_from_filename(fsa.file_name)
    if ctrl not in {"PK", "PK1", "PK2"}:
        return []

    assay = normalize_assay_qc(entry.get("assay", "UNKNOWN"))

    # Ladder-kanal iht master: ROX ladder = DATA4, LIZ ladder = DATA105
    ladder_channel = "DATA4" if entry.get("ladder") == "ROX" else "DATA105"

    wS = rules.sample_peak_window_bp
    wL = rules.ladder_peak_window_bp
    tuning = _merge_marker_tuning(wS, NON_SL_SAMPLE_TUNING.get(assay))
    sample_window_bp = float(tuning["sample_window_bp"])
    fallback_window_bp = float(tuning["fallback_window_bp"])

    # ---------------- SL (MNC_100 / MNC_20 etc): track bins 100/200/300/400/600 ----------------
    if assay == "SL":
        return [
            {"name": "SL_100", "kind": "sample", "expected_bp": 100.0, "channel": "primary", "window_bp": 20.0},
            {"name": "SL_200", "kind": "sample", "expected_bp": 200.0, "channel": "primary", "window_bp": 20.0},
            {"name": "SL_300", "kind": "sample", "expected_bp": 300.0, "channel": "primary", "window_bp": 20.0},
            {"name": "SL_400", "kind": "sample", "expected_bp": 400.0, "channel": "primary", "window_bp": 20.0},
            {"name": "SL_600", "kind": "sample", "expected_bp": 600.0, "channel": "primary", "window_bp": 40.0},
        ]

    # ---------------- FR1/FR2 (B=DATA1) + ladder ved 280 ----------------
    if assay == "FR1":
        return [
            {"name": "FR1_PK_DATA1_325", "kind": "sample", "expected_bp": 325.0, "channel": "DATA1", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "FR2":
        return [
            {"name": "FR2_PK_DATA1_260", "kind": "sample", "expected_bp": 260.0, "channel": "DATA1", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- FR3 (G=DATA2) + ladder ved 280 ----------------
    if assay == "FR3":
        return [
            {"name": "FR3_PK_DATA2_145", "kind": "sample", "expected_bp": 145.0, "channel": "DATA2", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- DHJH_D (G=DATA2) + ladder ved 300 ----------------
    if assay == "DHJH_D":
        return [
            {"name": "DHJH_D_PK_DATA2_139", "kind": "sample", "expected_bp": 139.0, "channel": "DATA2", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_300", "kind": "ladder", "expected_bp": 300.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- DHJH_E (B=DATA1) + ladder ved 150 ----------------
    if assay == "DHJH_E":
        return [
            {"name": "DHJH_E_PK_DATA1_109", "kind": "sample", "expected_bp": 109.0, "channel": "DATA1", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_150", "kind": "ladder", "expected_bp": 150.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- IGK: 279 (B=DATA1) og 150 (G=DATA2) + ladder ved 200 ----------------
    if assay == "IGK":
        return [
            {"name": "IGK_PK_DATA1_279", "kind": "sample", "expected_bp": 279.0, "channel": "DATA1", "window_bp": sample_window_bp},
            {"name": "IGK_PK_DATA2_150", "kind": "sample", "expected_bp": 150.0, "channel": "DATA2", "window_bp": sample_window_bp},
            {"name": "LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- KDE: (DATA3 i master) + ladder ved 200 ---------------- [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if assay == "KDE":
        return [
            {"name": "KDE_PK_DATA3_287", "kind": "sample", "expected_bp": 287.0, "channel": "DATA3", "window_bp": sample_window_bp},
            {"name": "KDE_PK_DATA3_377", "kind": "sample", "expected_bp": 377.0, "channel": "DATA3", "window_bp": sample_window_bp},
            {"name": "LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- TCRb: A/B/C (multi-kanal DATA1+DATA2) + ladder ved 280 ---------------- [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if assay == "TCRbA":
        return [
            {"name": "TCRbA_PK_DATA2_265", "kind": "sample", "expected_bp": 265.0, "channel": "DATA2", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "TCRbB":
        return [
            {"name": "TCRbB_PK_DATA1_254", "kind": "sample", "expected_bp": 254.0, "channel": "DATA1", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "TCRbC":
        return [
            {"name": "TCRbC_PK_DATA2_311", "kind": "sample", "expected_bp": 311.0, "channel": "DATA2", "window_bp": sample_window_bp},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
# ---------------- TCRg: A/B (multi-kanal DATA1+DATA2) + ladder ved 200 ----------------
    if assay == "TCRgA":
        if ctrl == "PK1":
            return [
                # PK1: B(blå)=DATA1 249, G(grønn)=DATA2 212, O*=ladder 200
                {"name": "TCRgA_PK1_DATA1_249", "kind": "sample", "expected_bp": 249.0, "channel": "DATA1", "window_bp": sample_window_bp},
                {"name": "TCRgA_PK1_DATA2_212", "kind": "sample", "expected_bp": 212.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgA_PK1_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        elif ctrl == "PK2":
            return [
                # PK2: G(grønn)=DATA2 163, O*=ladder 200
                {"name": "TCRgA_PK2_DATA2_163", "kind": "sample", "expected_bp": 163.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgA_PK2_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        else:
            # Hvis noen kjører "PK_TCRg_A" uten 1/2 (sjeldent), fall back:
            return [
                {"name": "TCRgA_PK_DATA1_249", "kind": "sample", "expected_bp": 249.0, "channel": "DATA1", "window_bp": sample_window_bp},
                {"name": "TCRgA_PK_DATA2_212", "kind": "sample", "expected_bp": 212.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgA_PK_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]

    if assay == "TCRgB":
        if ctrl == "PK1":
            return [
                # PK1: G(grønn)=DATA2 115, O*=ladder 200
                {"name": "TCRgB_PK1_DATA2_115", "kind": "sample", "expected_bp": 115.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgB_PK1_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        elif ctrl == "PK2":
            return [
                # PK2: G(grønn)=DATA2 178, O*=ladder 200
                {"name": "TCRgB_PK2_DATA2_178", "kind": "sample", "expected_bp": 178.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgB_PK2_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        else:
            return [
                {"name": "TCRgB_PK_DATA2_115", "kind": "sample", "expected_bp": 115.0, "channel": "DATA2", "window_bp": sample_window_bp},
                {"name": "TCRgB_PK_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]

    return []
# ======================================================================
# Peak finding nær forventet bp (height + area + found_bp)
# ======================================================================
def find_peak_near_bp(
    fsa,
    channel: str,
    target_bp: float,
    window_bp: float,
    baseline_correct: bool = True,
    name: str | None = None,
):
    """Finn maks i bp-vindu rundt target_bp i gitt kanal."""
    return _direct_peak_candidate_near_bp(
        fsa=fsa,
        channel=channel,
        target_bp=target_bp,
        window_bp=window_bp,
        baseline_correct=baseline_correct,
    )


def _direct_peak_candidate_near_bp(
    fsa,
    channel: str,
    target_bp: float,
    window_bp: float,
    baseline_correct: bool = True,
):
    """Returner enkeltkandidat for sterkeste peak i vinduet."""
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return {"ok": False, "reason": "mangler sample_data_with_basepairs"}

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return {"ok": False, "reason": "mangler time/basepairs kolonner"}

    if channel not in fsa.fsa:
        # Fallback for weird channel names or missing channels
        return {"ok": False, "reason": f"Channel {channel} not found in FSA file"}

    bp = raw_df["basepairs"].to_numpy()
    t = raw_df["time"].astype(int).to_numpy()

    win = (bp >= (target_bp - window_bp)) & (bp <= (target_bp + window_bp))
    if not np.any(win):
        return {"ok": False, "reason": "ingen punkter i bp-vindu"}

    bpw = bp[win]
    tw = t[win]
    trace = np.asarray(fsa.fsa[channel])

    valid = (tw >= 0) & (tw < len(trace))
    if not np.any(valid):
        return {"ok": False, "reason": "ingen gyldige time-indekser"}

    bpw = bpw[valid]
    tw = tw[valid]
    y = trace[tw].astype(float)

    # (valgfritt) samme baseline-korreksjon som master bruker i Plotly-plottene. [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if baseline_correct:
        try:
            baseline = estimate_running_baseline(trace, bin_size=200, quantile=0.10)
            y_full = trace.astype(float) - baseline
            y_full[y_full < 0] = 0.0
            y = y_full[tw]
        except Exception:
            pass

    if y.size == 0 or not np.any(np.isfinite(y)):
        return {"ok": False, "reason": "tom/NaN i vindu"}

    j = int(np.nanargmax(y))
    found_bp = float(bpw[j])
    height = float(y[j])
    area = float(np.nansum(y))

    return {
        "ok": True,
        "found_bp": found_bp,
        "height": height,
        "area": area,
        "search_window_bp": float(window_bp),
        "search_mode": "primary",
    }


def find_local_peak_candidates_near_bp(
    fsa,
    channel: str,
    target_bp: float,
    window_bp: float,
    baseline_correct: bool = True,
):
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return [{"ok": False, "reason": "mangler sample_data_with_basepairs"}]

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return [{"ok": False, "reason": "mangler time/basepairs kolonner"}]

    if channel not in fsa.fsa:
        return [{"ok": False, "reason": f"Channel {channel} not found in FSA file"}]

    bp = raw_df["basepairs"].to_numpy()
    t = raw_df["time"].astype(int).to_numpy()
    win = (bp >= (target_bp - window_bp)) & (bp <= (target_bp + window_bp))
    if not np.any(win):
        return [{"ok": False, "reason": "ingen punkter i bp-vindu"}]

    bpw = bp[win]
    tw = t[win]
    trace = np.asarray(fsa.fsa[channel])
    valid = (tw >= 0) & (tw < len(trace))
    if not np.any(valid):
        return [{"ok": False, "reason": "ingen gyldige time-indekser"}]

    bpw = bpw[valid]
    tw = tw[valid]
    y = trace[tw].astype(float)
    if baseline_correct:
        try:
            baseline = estimate_running_baseline(trace, bin_size=200, quantile=0.10)
            y_full = trace.astype(float) - baseline
            y_full[y_full < 0] = 0.0
            y = y_full[tw]
        except Exception:
            pass

    if y.size == 0 or not np.any(np.isfinite(y)):
        return [{"ok": False, "reason": "tom/NaN i vindu"}]

    area = float(np.nansum(y))
    candidates: list[dict] = []
    for idx in range(y.size):
        left = y[idx - 1] if idx > 0 else float("-inf")
        right = y[idx + 1] if idx + 1 < y.size else float("-inf")
        if not (y[idx] >= left and y[idx] >= right):
            continue
        candidates.append(
            {
                "ok": True,
                "found_bp": float(bpw[idx]),
                "height": float(y[idx]),
                "area": area,
                "search_window_bp": float(window_bp),
                "search_mode": "primary",
            }
        )

    if candidates:
        return candidates

    j = int(np.nanargmax(y))
    return [
        {
            "ok": True,
            "found_bp": float(bpw[j]),
            "height": float(y[j]),
            "area": area,
            "search_window_bp": float(window_bp),
            "search_mode": "primary",
        }
    ]


def _select_best_peak_candidate(
    candidates: list[dict],
    *,
    target_bp: float,
    reference_window_bp: float,
    scorer,
) -> dict:
    ok_candidates = [dict(candidate) for candidate in candidates if candidate.get("ok", False)]
    if not ok_candidates:
        return dict(candidates[0]) if candidates else {"ok": False, "reason": "ingen kandidater"}

    reference_height = max(float(candidate.get("height", 0.0) or 0.0) for candidate in ok_candidates)
    best = None
    best_score = float("-inf")
    for candidate in ok_candidates:
        score = scorer(candidate, target_bp, reference_window_bp, reference_height)
        candidate["selection_score"] = score
        if score > best_score:
            best = candidate
            best_score = score
    return best or ok_candidates[0]


def find_peak_near_bp_with_fallback(
    fsa,
    channel: str,
    target_bp: float,
    window_bp: float,
    fallback_window_bp: float | None = None,
    baseline_correct: bool = True,
    name: str | None = None,
):
    """Retry with a wider window when the peak is missed or sits on the window edge."""
    evaluation = evaluate_peak_near_bp_with_fallback(
        fsa=fsa,
        channel=channel,
        target_bp=target_bp,
        window_bp=window_bp,
        fallback_window_bp=fallback_window_bp,
        baseline_correct=baseline_correct,
        name=name,
    )
    return dict(evaluation["selected"])


def evaluate_peak_near_bp_with_fallback(
    fsa,
    channel: str,
    target_bp: float,
    window_bp: float,
    fallback_window_bp: float | None = None,
    baseline_correct: bool = True,
    name: str | None = None,
):
    """Return both evaluated attempts and the selected result for PK tracking."""
    tuning = _tuning_for_marker_target(target_bp, window_bp)
    primary_window_bp = max(float(window_bp), float(tuning["sample_window_bp"]))
    fallback_window_bp = float(fallback_window_bp or window_bp)
    fallback_window_bp = max(fallback_window_bp, float(tuning["fallback_window_bp"]), primary_window_bp)
    selection_mode = str(tuning.get("selection_mode", "legacy"))

    if selection_mode == "local_score":
        primary = _select_best_peak_candidate(
            find_local_peak_candidates_near_bp(
                fsa=fsa,
                channel=channel,
                target_bp=target_bp,
                window_bp=primary_window_bp,
                baseline_correct=baseline_correct,
            ),
            target_bp=target_bp,
            reference_window_bp=fallback_window_bp,
            scorer=_local_candidate_score,
        )
    else:
        primary = _direct_peak_candidate_near_bp(
            fsa=fsa,
            channel=channel,
            target_bp=target_bp,
            window_bp=primary_window_bp,
            baseline_correct=baseline_correct,
        )
    primary = dict(primary)
    primary.setdefault("search_mode", "primary")
    if primary.get("search_window_bp") is None:
        primary["search_window_bp"] = float(primary_window_bp)
    primary["delta_bp"] = (
        float(primary["found_bp"]) - float(target_bp)
        if primary.get("ok", False) and primary.get("found_bp") is not None
        else None
    )

    candidates = [primary]
    if fallback_window_bp <= float(primary_window_bp):
        return {
            "selected": primary,
            "selected_index": 0,
            "candidates": candidates,
            "tuning": {
                "selection_mode": selection_mode,
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    if selection_mode == "local_score":
        fallback = _select_best_peak_candidate(
            find_local_peak_candidates_near_bp(
                fsa=fsa,
                channel=channel,
                target_bp=target_bp,
                window_bp=fallback_window_bp,
                baseline_correct=baseline_correct,
            ),
            target_bp=target_bp,
            reference_window_bp=fallback_window_bp,
            scorer=_local_candidate_score,
        )
    else:
        fallback = _direct_peak_candidate_near_bp(
            fsa=fsa,
            channel=channel,
            target_bp=target_bp,
            window_bp=fallback_window_bp,
            baseline_correct=baseline_correct,
        )
    fallback = dict(fallback)
    fallback["search_mode"] = "fallback"
    fallback["search_window_bp"] = float(fallback_window_bp)
    fallback["delta_bp"] = (
        float(fallback["found_bp"]) - float(target_bp)
        if fallback.get("ok", False) and fallback.get("found_bp") is not None
        else None
    )
    fallback["fallback_from_window_bp"] = float(primary_window_bp)
    
    candidates.append(fallback)
    if not fallback.get("ok", False):
        return {
            "selected": primary,
            "selected_index": 0,
            "candidates": candidates,
            "tuning": {
                "selection_mode": selection_mode,
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    if not primary.get("ok", False):
        fallback["fallback_from_window_bp"] = float(primary_window_bp)
        return {
            "selected": fallback,
            "selected_index": 1,
            "candidates": candidates,
            "tuning": {
                "selection_mode": selection_mode,
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    if selection_mode in {"score", "local_score"}:
        reference_height = max(
            float(primary.get("height", 0.0) or 0.0),
            float(fallback.get("height", 0.0) or 0.0),
            1.0,
        )
        scorer = _local_candidate_score if selection_mode == "local_score" else _candidate_score
        primary_score = scorer(primary, target_bp, fallback_window_bp, reference_height)
        fallback_score = scorer(fallback, target_bp, fallback_window_bp, reference_height)
        primary["selection_score"] = primary_score
        fallback["selection_score"] = fallback_score
        if fallback_score > (primary_score + 0.05):
            fallback["fallback_from_window_bp"] = float(primary_window_bp)
            return {
                "selected": fallback,
                "selected_index": 1,
                "candidates": candidates,
                "tuning": {
                    "selection_mode": selection_mode,
                    "primary_window_bp": float(primary_window_bp),
                    "fallback_window_bp": float(fallback_window_bp),
                },
            }
        return {
            "selected": primary,
            "selected_index": 0,
            "candidates": candidates,
            "tuning": {
                "selection_mode": selection_mode,
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    primary_delta = abs(float(primary["found_bp"]) - float(target_bp))
    fallback_delta = abs(float(fallback["found_bp"]) - float(target_bp))
    if fallback_delta < primary_delta:
        fallback["fallback_from_window_bp"] = float(primary_window_bp)
        return {
            "selected": fallback,
            "selected_index": 1,
            "candidates": candidates,
            "tuning": {
                "selection_mode": "legacy",
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    primary_height = float(primary.get("height", 0.0) or 0.0)
    fallback_height = float(fallback.get("height", 0.0) or 0.0)
    if primary_height <= 0:
        primary_height = 1.0
    primary_near_edge = primary_delta >= max(float(primary_window_bp) * 0.8, float(primary_window_bp) - 0.5)
    fallback_not_much_further = fallback_delta <= (primary_delta + 1.0)
    if fallback_height >= (primary_height * 1.5) and (primary_near_edge or fallback_not_much_further):
        fallback["fallback_from_window_bp"] = float(primary_window_bp)
        fallback["search_mode"] = "fallback"
        return {
            "selected": fallback,
            "selected_index": 1,
            "candidates": candidates,
            "tuning": {
                "selection_mode": "legacy",
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    if fallback_height >= (primary_height * 2.0):
        fallback["fallback_from_window_bp"] = float(primary_window_bp)
        fallback["search_mode"] = "fallback"
        return {
            "selected": fallback,
            "selected_index": 1,
            "candidates": candidates,
            "tuning": {
                "selection_mode": "legacy",
                "primary_window_bp": float(primary_window_bp),
                "fallback_window_bp": float(fallback_window_bp),
            },
        }

    return {
        "selected": primary,
        "selected_index": 0,
        "candidates": candidates,
        "tuning": {
            "selection_mode": "legacy",
            "primary_window_bp": float(primary_window_bp),
            "fallback_window_bp": float(fallback_window_bp),
        },
    }
