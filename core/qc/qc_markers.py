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
CONTROL_PREFIX_RE = re.compile(r"^(PK1|PK2|PK|NK|RK)_", re.IGNORECASE) # PK1_TCRgA...

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


STAGE_PREFIX_RE = re.compile(r"^\d{5}_[a-f0-9]{8}_", re.IGNORECASE)

def strip_stage_prefix(name: str) -> str:
    return STAGE_PREFIX_RE.sub("", name)

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
            {"name": "FR1_PK_DATA1_325", "kind": "sample", "expected_bp": 325.0, "channel": "DATA1", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "FR2":
        return [
            {"name": "FR2_PK_DATA1_260", "kind": "sample", "expected_bp": 260.0, "channel": "DATA1", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- FR3 (G=DATA2) + ladder ved 280 ----------------
    if assay == "FR3":
        return [
            {"name": "FR3_PK_DATA2_145", "kind": "sample", "expected_bp": 145.0, "channel": "DATA2", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- DHJH_D (G=DATA2) + ladder ved 300 ----------------
    if assay == "DHJH_D":
        return [
            {"name": "DHJH_D_PK_DATA2_139", "kind": "sample", "expected_bp": 139.0, "channel": "DATA2", "window_bp": wS},
            {"name": "ROX_Ladder_300", "kind": "ladder", "expected_bp": 300.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- DHJH_E (B=DATA1) + ladder ved 150 ----------------
    if assay == "DHJH_E":
        return [
            {"name": "DHJH_E_PK_DATA1_109", "kind": "sample", "expected_bp": 109.0, "channel": "DATA1", "window_bp": wS},
            {"name": "ROX_Ladder_150", "kind": "ladder", "expected_bp": 150.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- IGK: 279 (B=DATA1) og 150 (G=DATA2) + ladder ved 200 ----------------
    if assay == "IGK":
        return [
            {"name": "IGK_PK_DATA1_279", "kind": "sample", "expected_bp": 279.0, "channel": "DATA1", "window_bp": wS},
            {"name": "IGK_PK_DATA2_150", "kind": "sample", "expected_bp": 150.0, "channel": "DATA2", "window_bp": wS},
            {"name": "LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- KDE: (DATA3 i master) + ladder ved 200 ---------------- [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if assay == "KDE":
        return [
            {"name": "KDE_PK_DATA3_287", "kind": "sample", "expected_bp": 287.0, "channel": "DATA3", "window_bp": wS},
            {"name": "KDE_PK_DATA3_377", "kind": "sample", "expected_bp": 377.0, "channel": "DATA3", "window_bp": wS},
            {"name": "LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
        ]

    # ---------------- TCRb: A/B/C (multi-kanal DATA1+DATA2) + ladder ved 280 ---------------- [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if assay == "TCRbA":
        return [
            {"name": "TCRbA_PK_DATA2_265", "kind": "sample", "expected_bp": 265.0, "channel": "DATA2", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "TCRbB":
        return [
            {"name": "TCRbB_PK_DATA1_254", "kind": "sample", "expected_bp": 254.0, "channel": "DATA1", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
    if assay == "TCRbC":
        return [
            {"name": "TCRbC_PK_DATA2_311", "kind": "sample", "expected_bp": 311.0, "channel": "DATA2", "window_bp": wS},
            {"name": "ROX_Ladder_280", "kind": "ladder", "expected_bp": 280.0, "channel": ladder_channel, "window_bp": wL},
        ]
# ---------------- TCRg: A/B (multi-kanal DATA1+DATA2) + ladder ved 200 ----------------
    if assay == "TCRgA":
        if ctrl == "PK1":
            return [
                # PK1: B(blå)=DATA1 249, G(grønn)=DATA2 212, O*=ladder 200
                {"name": "TCRgA_PK1_DATA1_249", "kind": "sample", "expected_bp": 249.0, "channel": "DATA1", "window_bp": wS},
                {"name": "TCRgA_PK1_DATA2_212", "kind": "sample", "expected_bp": 212.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgA_PK1_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        elif ctrl == "PK2":
            return [
                # PK2: G(grønn)=DATA2 163, O*=ladder 200
                {"name": "TCRgA_PK2_DATA2_163", "kind": "sample", "expected_bp": 163.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgA_PK2_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        else:
            # Hvis noen kjører "PK_TCRg_A" uten 1/2 (sjeldent), fall back:
            return [
                {"name": "TCRgA_PK_DATA1_249", "kind": "sample", "expected_bp": 249.0, "channel": "DATA1", "window_bp": wS},
                {"name": "TCRgA_PK_DATA2_212", "kind": "sample", "expected_bp": 212.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgA_PK_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]

    if assay == "TCRgB":
        if ctrl == "PK1":
            return [
                # PK1: G(grønn)=DATA2 115, O*=ladder 200
                {"name": "TCRgB_PK1_DATA2_115", "kind": "sample", "expected_bp": 115.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgB_PK1_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        elif ctrl == "PK2":
            return [
                # PK2: G(grønn)=DATA2 178, O*=ladder 200
                {"name": "TCRgB_PK2_DATA2_178", "kind": "sample", "expected_bp": 178.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgB_PK2_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]
        else:
            return [
                {"name": "TCRgB_PK_DATA2_115", "kind": "sample", "expected_bp": 115.0, "channel": "DATA2", "window_bp": wS},
                {"name": "TCRgB_PK_LIZ_Ladder_200", "kind": "ladder", "expected_bp": 200.0, "channel": ladder_channel, "window_bp": wL},
            ]

    return []
# ======================================================================
# Peak finding nær forventet bp (height + area + found_bp)
# ======================================================================
def find_peak_near_bp(fsa, channel: str, target_bp: float, window_bp: float, baseline_correct: bool = True):
    """
    Finn maks i bp-vindu rundt target_bp i gitt kanal.
    Returnerer dict med found_bp, height, area og ok-flag.
    """
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return {"ok": False, "reason": "mangler sample_data_with_basepairs"}

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return {"ok": False, "reason": "mangler time/basepairs kolonner"}

    if channel not in fsa.fsa:
        return {"ok": False, "reason": f"kanal {channel} finnes ikke i fsa"}

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

    return {"ok": True, "found_bp": found_bp, "height": height, "area": area}
