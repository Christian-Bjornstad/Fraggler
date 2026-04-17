from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.chart import BarChart, Reference
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter


HEADER_FILL = PatternFill("solid", fgColor="1F4E78")
CARD_BLUE = PatternFill("solid", fgColor="D9EAF7")
CARD_TEAL = PatternFill("solid", fgColor="DDF4F1")
CARD_ORANGE = PatternFill("solid", fgColor="FCE4D6")
CARD_RED = PatternFill("solid", fgColor="FBE5E7")
CARD_GOLD = PatternFill("solid", fgColor="FFF2CC")
HEADER_FONT = Font(color="FFFFFF", bold=True)
BOLD = Font(bold=True)
THIN_GRAY = Side(style="thin", color="D9E2F2")
BOX_BORDER = Border(left=THIN_GRAY, right=THIN_GRAY, top=THIN_GRAY, bottom=THIN_GRAY)
DASHBOARD_SHEETS = [
    "Dashboard",
    "Dashboard_Data",
    "Assay_Summary",
    "Run_Summary",
    "Control_Summary",
    "PK_Sample_Delta",
    "PK_Ladder_Delta",
]


def refresh_clonality_tracking_dashboard(
    excel_path: Path,
    *,
    dashboard_title: str = "Fraggler Klonalitet Dashboard",
) -> None:
    if not excel_path.exists():
        return

    with pd.ExcelFile(excel_path, engine="openpyxl") as xls:
        required = {"Patient_Runs", "Control_Runs", "PK_Peaks"}
        if not required.issubset(set(xls.sheet_names)):
            return

    wb = load_workbook(excel_path)
    _ensure_abs_delta_column(wb["PK_Peaks"])
    for name in DASHBOARD_SHEETS:
        if name in wb.sheetnames:
            del wb[name]

    dashboard = wb.create_sheet("Dashboard", 0)
    data_ws = wb.create_sheet("Dashboard_Data")
    data_ws.sheet_state = "hidden"
    assay_ws = wb.create_sheet("Assay_Summary")
    run_ws = wb.create_sheet("Run_Summary")
    control_ws = wb.create_sheet("Control_Summary")
    pk_sample_ws = wb.create_sheet("PK_Sample_Delta")
    pk_ladder_ws = wb.create_sheet("PK_Ladder_Delta")

    try:
        wb.calculation.calcMode = "auto"
        wb.calculation.fullCalcOnLoad = True
        wb.calculation.forceFullCalc = True
    except Exception:
        pass

    patient_cols = _col_map(wb["Patient_Runs"])
    control_cols = _col_map(wb["Control_Runs"])
    pk_cols = _col_map(wb["PK_Peaks"])

    assays = _sorted_unique_values(wb["Patient_Runs"], patient_cols.get("Assay"), wb["Control_Runs"], control_cols.get("Assay"))
    controls = _sorted_unique_values(wb["Control_Runs"], control_cols.get("Control"))
    review_pairs = _control_assay_pairs(wb["Control_Runs"], control_cols)

    _write_helper_lists(data_ws, assays, controls, review_pairs)
    _build_assay_summary(assay_ws, assays, patient_cols, control_cols)
    _build_run_summary(run_ws, assays, patient_cols)
    _build_control_summary(control_ws, review_pairs, control_cols)
    _build_pk_summary(pk_sample_ws, assays, pk_cols, kind="sample")
    _build_pk_summary(pk_ladder_ws, assays, pk_cols, kind="ladder")
    _build_dashboard(
        dashboard,
        dashboard_title=dashboard_title,
        assay_count=len(assays),
        patient_cols=patient_cols,
        control_cols=control_cols,
        pk_cols=pk_cols,
    )
    _add_dashboard_charts(dashboard, assay_ws, pk_sample_ws)

    wb.save(excel_path)


def _ensure_abs_delta_column(ws) -> None:
    headers = [str(cell.value).strip() if cell.value is not None else "" for cell in ws[1]]
    if "AbsDeltaBP" in headers:
        return

    delta_idx = headers.index("DeltaBP") + 1 if "DeltaBP" in headers else None
    next_col = ws.max_column + 1
    ws.cell(1, next_col, "AbsDeltaBP")
    for row_idx in range(2, ws.max_row + 1):
        delta_value = ws.cell(row_idx, delta_idx).value if delta_idx else None
        if delta_value in ("", None):
            continue
        try:
            ws.cell(row_idx, next_col, abs(float(delta_value)))
        except (TypeError, ValueError):
            continue


def _col_map(ws) -> dict[str, str]:
    result: dict[str, str] = {}
    for cell in ws[1]:
        value = cell.value
        if value is None:
            continue
        result[str(value)] = get_column_letter(cell.column)
    return result


def _range_ref(sheet_name: str, col: str, *, start_row: int = 2, end_row: int = 1048576) -> str:
    return f"{sheet_name}!${col}${start_row}:${col}${end_row}"


def _style_table(ws, header_row: int, data_end_row: int, start_col: int, end_col: int) -> None:
    for cell in ws[header_row]:
        if start_col <= cell.column <= end_col:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
            cell.alignment = Alignment(horizontal="center", vertical="center")
            cell.border = BOX_BORDER
    for row in ws.iter_rows(min_row=header_row + 1, max_row=max(data_end_row, header_row + 1), min_col=start_col, max_col=end_col):
        for cell in row:
            cell.border = BOX_BORDER
            if cell.row % 2 == 0:
                cell.fill = PatternFill("solid", fgColor="F8FBFF")


def _autofit(ws) -> None:
    for column_cells in ws.columns:
        length = 0
        col_idx = column_cells[0].column
        for cell in column_cells:
            value = "" if cell.value is None else str(cell.value)
            length = max(length, len(value))
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max(length + 2, 10), 30)


def _sorted_unique_values(*pairs) -> list[str]:
    values: set[str] = set()
    for idx in range(0, len(pairs), 2):
        ws = pairs[idx]
        col = pairs[idx + 1] if idx + 1 < len(pairs) else None
        if ws is None or not col:
            continue
        for cell in ws[col][1:]:
            value = str(cell.value or "").strip()
            if value:
                values.add(value)
    return sorted(values)


def _control_assay_pairs(ws, cols: dict[str, str]) -> list[tuple[str, str]]:
    control_col = cols.get("Control")
    assay_col = cols.get("Assay")
    if not control_col or not assay_col:
        return []
    pairs: set[tuple[str, str]] = set()
    for row_idx in range(2, ws.max_row + 1):
        control = str(ws[f"{control_col}{row_idx}"].value or "").strip()
        assay = str(ws[f"{assay_col}{row_idx}"].value or "").strip()
        if control and assay:
            pairs.add((control, assay))
    return sorted(pairs, key=lambda item: (item[0], item[1]))


def _write_helper_lists(ws, assays: list[str], controls: list[str], review_pairs: list[tuple[str, str]]) -> None:
    ws["A1"] = "Assay"
    ws["E1"] = "Control"
    ws["F1"] = "Assay"
    for row_idx, assay in enumerate(assays, start=2):
        ws.cell(row_idx, 1, assay)
    for row_idx, control in enumerate(controls, start=2):
        ws.cell(row_idx, 5, control)
    for row_idx, (control, assay) in enumerate(review_pairs, start=2):
        ws.cell(row_idx, 5, control)
        ws.cell(row_idx, 6, assay)


def _build_assay_summary(ws, assays: list[str], patient_cols: dict[str, str], control_cols: dict[str, str]) -> None:
    headers = ["Assay", "Files", "PatientFiles", "ControlFiles", "LadderReview", "AvgR2", "PartialFits", "ReviewRate"]
    ws.append(headers)
    p_assay = patient_cols["Assay"]
    p_qc = patient_cols["LadderQC"]
    p_r2 = patient_cols["LadderR2"]
    p_expected = patient_cols["LadderExpectedStepCount"]
    p_fitted = patient_cols["LadderFittedStepCount"]
    c_assay = control_cols["Assay"]
    c_qc = control_cols["LadderQC"]
    c_r2 = control_cols["LadderR2"]
    c_expected = control_cols["LadderExpectedStepCount"]
    c_fitted = control_cols["LadderFittedStepCount"]
    for row_idx, assay in enumerate(assays, start=2):
        ws.cell(row_idx, 1, assay)
        ws.cell(row_idx, 2, f'=COUNTIF(Patient_Runs!${p_assay}:${p_assay},$A{row_idx})+COUNTIF(Control_Runs!${c_assay}:${c_assay},$A{row_idx})')
        ws.cell(row_idx, 3, f'=COUNTIF(Patient_Runs!${p_assay}:${p_assay},$A{row_idx})')
        ws.cell(row_idx, 4, f'=COUNTIF(Control_Runs!${c_assay}:${c_assay},$A{row_idx})')
        ws.cell(row_idx, 5, f'=COUNTIFS(Patient_Runs!${p_assay}:${p_assay},$A{row_idx},Patient_Runs!${p_qc}:${p_qc},"<>",Patient_Runs!${p_qc}:${p_qc},"<>ok")+COUNTIFS(Control_Runs!${c_assay}:${c_assay},$A{row_idx},Control_Runs!${c_qc}:${c_qc},"<>",Control_Runs!${c_qc}:${c_qc},"<>ok")')
        ws.cell(row_idx, 6, f'=IFERROR((SUMIF(Patient_Runs!${p_assay}:${p_assay},$A{row_idx},Patient_Runs!${p_r2}:${p_r2})+SUMIF(Control_Runs!${c_assay}:${c_assay},$A{row_idx},Control_Runs!${c_r2}:${c_r2}))/$B{row_idx},0)')
        ws.cell(row_idx, 7, f'=SUMPRODUCT(--(Patient_Runs!${p_assay}$2:${p_assay}$1048576=$A{row_idx}),--(Patient_Runs!${p_fitted}$2:${p_fitted}$1048576<Patient_Runs!${p_expected}$2:${p_expected}$1048576))+SUMPRODUCT(--(Control_Runs!${c_assay}$2:${c_assay}$1048576=$A{row_idx}),--(Control_Runs!${c_fitted}$2:${c_fitted}$1048576<Control_Runs!${c_expected}$2:${c_expected}$1048576))')
        ws.cell(row_idx, 8, f'=IFERROR(E{row_idx}/B{row_idx},0)')
    _style_table(ws, 1, max(len(assays) + 1, 2), 1, len(headers))
    ws.freeze_panes = "A2"
    _autofit(ws)


def _build_run_summary(ws, assays: list[str], patient_cols: dict[str, str]) -> None:
    headers = ["Assay", "Files", "ReviewFiles", "AvgR2", "PartialFits"]
    ws.append(headers)
    p_assay = patient_cols["Assay"]
    p_qc = patient_cols["LadderQC"]
    p_r2 = patient_cols["LadderR2"]
    p_expected = patient_cols["LadderExpectedStepCount"]
    p_fitted = patient_cols["LadderFittedStepCount"]
    for row_idx, assay in enumerate(assays, start=2):
        ws.cell(row_idx, 1, assay)
        ws.cell(row_idx, 2, f'=COUNTIF(Patient_Runs!${p_assay}:${p_assay},$A{row_idx})')
        ws.cell(row_idx, 3, f'=COUNTIFS(Patient_Runs!${p_assay}:${p_assay},$A{row_idx},Patient_Runs!${p_qc}:${p_qc},"<>",Patient_Runs!${p_qc}:${p_qc},"<>ok")')
        ws.cell(row_idx, 4, f'=IFERROR(SUMIF(Patient_Runs!${p_assay}:${p_assay},$A{row_idx},Patient_Runs!${p_r2}:${p_r2})/B{row_idx},0)')
        ws.cell(row_idx, 5, f'=SUMPRODUCT(--(Patient_Runs!${p_assay}$2:${p_assay}$1048576=$A{row_idx}),--(Patient_Runs!${p_fitted}$2:${p_fitted}$1048576<Patient_Runs!${p_expected}$2:${p_expected}$1048576))')
    _style_table(ws, 1, max(len(assays) + 1, 2), 1, len(headers))
    ws.freeze_panes = "A2"
    _autofit(ws)


def _build_control_summary(ws, pairs: list[tuple[str, str]], control_cols: dict[str, str]) -> None:
    headers = ["Control", "Assay", "Files", "ReviewFiles", "AvgR2", "PartialFits"]
    ws.append(headers)
    c_control = control_cols["Control"]
    c_assay = control_cols["Assay"]
    c_qc = control_cols["LadderQC"]
    c_r2 = control_cols["LadderR2"]
    c_expected = control_cols["LadderExpectedStepCount"]
    c_fitted = control_cols["LadderFittedStepCount"]
    for row_idx, (control, assay) in enumerate(pairs, start=2):
        ws.cell(row_idx, 1, control)
        ws.cell(row_idx, 2, assay)
        ws.cell(row_idx, 3, f'=COUNTIFS(Control_Runs!${c_control}:${c_control},$A{row_idx},Control_Runs!${c_assay}:${c_assay},$B{row_idx})')
        ws.cell(row_idx, 4, f'=COUNTIFS(Control_Runs!${c_control}:${c_control},$A{row_idx},Control_Runs!${c_assay}:${c_assay},$B{row_idx},Control_Runs!${c_qc}:${c_qc},"<>",Control_Runs!${c_qc}:${c_qc},"<>ok")')
        ws.cell(row_idx, 5, f'=IFERROR(SUMIFS(Control_Runs!${c_r2}:${c_r2},Control_Runs!${c_control}:${c_control},$A{row_idx},Control_Runs!${c_assay}:${c_assay},$B{row_idx})/C{row_idx},0)')
        ws.cell(row_idx, 6, f'=SUMPRODUCT(--(Control_Runs!${c_control}$2:${c_control}$1048576=$A{row_idx}),--(Control_Runs!${c_assay}$2:${c_assay}$1048576=$B{row_idx}),--(Control_Runs!${c_fitted}$2:${c_fitted}$1048576<Control_Runs!${c_expected}$2:${c_expected}$1048576))')
    _style_table(ws, 1, max(len(pairs) + 1, 2), 1, len(headers))
    ws.freeze_panes = "A2"
    _autofit(ws)


def _build_pk_summary(ws, assays: list[str], pk_cols: dict[str, str], *, kind: str) -> None:
    headers = ["Assay", "MarkerRows", "MeanAbsDeltaBP", "MaxAbsDeltaBP", "Over2bp", "Over5bp", "AvgHeight"]
    ws.append(headers)
    pk_assay = pk_cols["Assay"]
    pk_kind = pk_cols["Kind"]
    pk_abs = pk_cols["AbsDeltaBP"]
    pk_height = pk_cols["Height"]
    for row_idx, assay in enumerate(assays, start=2):
        ws.cell(row_idx, 1, assay)
        ws.cell(row_idx, 2, f'=COUNTIFS(PK_Peaks!${pk_assay}:${pk_assay},$A{row_idx},PK_Peaks!${pk_kind}:${pk_kind},"{kind}")')
        ws.cell(row_idx, 3, f'=IFERROR(AVERAGEIFS(PK_Peaks!${pk_abs}:${pk_abs},PK_Peaks!${pk_assay}:${pk_assay},$A{row_idx},PK_Peaks!${pk_kind}:${pk_kind},"{kind}"),0)')
        ws.cell(row_idx, 4, f'=IFERROR(MAXIFS(PK_Peaks!${pk_abs}:${pk_abs},PK_Peaks!${pk_assay}:${pk_assay},$A{row_idx},PK_Peaks!${pk_kind}:${pk_kind},"{kind}"),0)')
        ws.cell(row_idx, 5, f'=SUMPRODUCT(--(PK_Peaks!${pk_assay}$2:${pk_assay}$1048576=$A{row_idx}),--(PK_Peaks!${pk_kind}$2:${pk_kind}$1048576="{kind}"),--(PK_Peaks!${pk_abs}$2:${pk_abs}$1048576>2))')
        ws.cell(row_idx, 6, f'=SUMPRODUCT(--(PK_Peaks!${pk_assay}$2:${pk_assay}$1048576=$A{row_idx}),--(PK_Peaks!${pk_kind}$2:${pk_kind}$1048576="{kind}"),--(PK_Peaks!${pk_abs}$2:${pk_abs}$1048576>5))')
        ws.cell(row_idx, 7, f'=IFERROR(AVERAGEIFS(PK_Peaks!${pk_height}:${pk_height},PK_Peaks!${pk_assay}:${pk_assay},$A{row_idx},PK_Peaks!${pk_kind}:${pk_kind},"{kind}"),0)')
    _style_table(ws, 1, max(len(assays) + 1, 2), 1, len(headers))
    ws.freeze_panes = "A2"
    _autofit(ws)


def _build_dashboard(ws, *, dashboard_title: str, assay_count: int, patient_cols: dict[str, str], control_cols: dict[str, str], pk_cols: dict[str, str]) -> None:
    ws.sheet_view.showGridLines = False
    for col, width in {"A": 18, "B": 16, "C": 16, "D": 16, "E": 18, "F": 18, "G": 18, "H": 18, "I": 4, "J": 16, "K": 16, "L": 16, "M": 16, "N": 16, "O": 16, "P": 16, "Q": 16}.items():
        ws.column_dimensions[col].width = width

    ws.merge_cells("A1:Q2")
    ws["A1"] = dashboard_title
    ws["A1"].font = Font(size=18, bold=True, color="FFFFFF")
    ws["A1"].fill = HEADER_FILL
    ws["A1"].alignment = Alignment(horizontal="center", vertical="center")

    labels = ["Tracked Runs", "Unique Assays", "Review Files", "Ladder OK Rate", "PK Marker Rows", "PK Sample Mean |delta|", "PK Sample >2 bp", "PK Ladder Mean |delta|"]
    fills = [CARD_BLUE, CARD_TEAL, CARD_RED, CARD_GOLD, CARD_BLUE, CARD_ORANGE, CARD_RED, CARD_TEAL]
    for idx, (label, fill) in enumerate(zip(labels, fills), start=1):
        cell = ws.cell(4, idx, label)
        cell.fill = fill
        cell.font = BOLD
        cell.alignment = Alignment(horizontal="center")
        cell.border = BOX_BORDER
    p_id = patient_cols["IdentityKey"]
    p_qc = patient_cols["LadderQC"]
    p_r2 = patient_cols["LadderR2"]
    p_expected = patient_cols["LadderExpectedStepCount"]
    p_fitted = patient_cols["LadderFittedStepCount"]
    c_id = control_cols["IdentityKey"]
    c_qc = control_cols["LadderQC"]
    c_r2 = control_cols["LadderR2"]
    c_expected = control_cols["LadderExpectedStepCount"]
    c_fitted = control_cols["LadderFittedStepCount"]
    pk_id = pk_cols["IdentityKey"]
    pk_kind = pk_cols["Kind"]
    pk_abs = pk_cols["AbsDeltaBP"]
    pk_assay = pk_cols["Assay"]
    ws["A5"] = f'=COUNTA(Patient_Runs!${p_id}:${p_id})+COUNTA(Control_Runs!${c_id}:${c_id})-2'
    ws["B5"] = f"={assay_count}"
    ws["C5"] = f'=COUNTIFS(Patient_Runs!${p_qc}:${p_qc},"<>",Patient_Runs!${p_qc}:${p_qc},"<>ok")+COUNTIFS(Control_Runs!${c_qc}:${c_qc},"<>",Control_Runs!${c_qc}:${c_qc},"<>ok")'
    ws["D5"] = f'=IFERROR((COUNTIF(Patient_Runs!${p_qc}:${p_qc},"ok")+COUNTIF(Control_Runs!${c_qc}:${c_qc},"ok"))/(COUNTA(Patient_Runs!${p_id}:${p_id})+COUNTA(Control_Runs!${c_id}:${c_id})-2),0)'
    ws["E5"] = f'=COUNTA(PK_Peaks!${pk_id}:${pk_id})-1'
    ws["F5"] = f'=IFERROR(AVERAGEIFS(PK_Peaks!${pk_abs}:${pk_abs},PK_Peaks!${pk_kind}:${pk_kind},"sample",PK_Peaks!${pk_assay}:${pk_assay},"<>SL"),0)'
    ws["G5"] = f'=COUNTIFS(PK_Peaks!${pk_kind}:${pk_kind},"sample",PK_Peaks!${pk_abs}:${pk_abs},">2",PK_Peaks!${pk_assay}:${pk_assay},"<>SL")'
    ws["H5"] = f'=IFERROR(AVERAGEIFS(PK_Peaks!${pk_abs}:${pk_abs},PK_Peaks!${pk_kind}:${pk_kind},"ladder"),0)'
    for col_idx in range(1, 9):
        cell = ws.cell(5, col_idx)
        cell.border = BOX_BORDER
        cell.alignment = Alignment(horizontal="center")

    ws["A7"] = "Ladder Overview"
    ws["A7"].font = Font(size=13, bold=True)
    ws["J7"] = "PK Sample Delta Focus"
    ws["J7"].font = Font(size=13, bold=True)
    overview_headers = ["Scope", "Runs", "Ladder OK", "Review Required", "OK Rate", "Avg R2", "Partial Fits"]
    for idx, header in enumerate(overview_headers, start=1):
        ws.cell(8, idx, header)
    ws["B9"] = f'=COUNTA(Patient_Runs!${p_id}:${p_id})+COUNTA(Control_Runs!${c_id}:${c_id})-2'
    ws["C9"] = f'=COUNTIF(Patient_Runs!${p_qc}:${p_qc},"ok")+COUNTIF(Control_Runs!${c_qc}:${c_qc},"ok")'
    ws["D9"] = f'=COUNTIFS(Patient_Runs!${p_qc}:${p_qc},"<>",Patient_Runs!${p_qc}:${p_qc},"<>ok")+COUNTIFS(Control_Runs!${c_qc}:${c_qc},"<>",Control_Runs!${c_qc}:${c_qc},"<>ok")'
    ws["E9"] = '=IFERROR(C9/B9,0)'
    ws["F9"] = f'=IFERROR((SUM(Patient_Runs!${p_r2}:${p_r2})+SUM(Control_Runs!${c_r2}:${c_r2}))/B9,0)'
    ws["G9"] = f'=SUMPRODUCT(--(Patient_Runs!${p_fitted}$2:${p_fitted}$1048576<Patient_Runs!${p_expected}$2:${p_expected}$1048576))+SUMPRODUCT(--(Control_Runs!${c_fitted}$2:${c_fitted}$1048576<Control_Runs!${c_expected}$2:${c_expected}$1048576))'
    ws["A9"] = "All"
    ws["A10"] = "Patient"
    ws["B10"] = f'=COUNTA(Patient_Runs!${p_id}:${p_id})-1'
    ws["C10"] = f'=COUNTIF(Patient_Runs!${p_qc}:${p_qc},"ok")'
    ws["D10"] = f'=COUNTIFS(Patient_Runs!${p_qc}:${p_qc},"<>",Patient_Runs!${p_qc}:${p_qc},"<>ok")'
    ws["E10"] = '=IFERROR(C10/B10,0)'
    ws["F10"] = f'=IFERROR(SUM(Patient_Runs!${p_r2}:${p_r2})/B10,0)'
    ws["G10"] = f'=SUMPRODUCT(--(Patient_Runs!${p_fitted}$2:${p_fitted}$1048576<Patient_Runs!${p_expected}$2:${p_expected}$1048576))'
    ws["A11"] = "Control"
    ws["B11"] = f'=COUNTA(Control_Runs!${c_id}:${c_id})-1'
    ws["C11"] = f'=COUNTIF(Control_Runs!${c_qc}:${c_qc},"ok")'
    ws["D11"] = f'=COUNTIFS(Control_Runs!${c_qc}:${c_qc},"<>",Control_Runs!${c_qc}:${c_qc},"<>ok")'
    ws["E11"] = '=IFERROR(C11/B11,0)'
    ws["F11"] = f'=IFERROR(SUM(Control_Runs!${c_r2}:${c_r2})/B11,0)'
    ws["G11"] = f'=SUMPRODUCT(--(Control_Runs!${c_fitted}$2:${c_fitted}$1048576<Control_Runs!${c_expected}$2:${c_expected}$1048576))'
    _style_table(ws, 8, 11, 1, 7)

    focus_headers = ["Assay", "MarkerRows", "MeanAbsDeltaBP", "MaxAbsDeltaBP", "Over2bp", "Over5bp", "AvgHeight"]
    for idx, header in enumerate(focus_headers, start=10):
        ws.cell(8, idx, header)
    for offset in range(0, min(assay_count, 12)):
        src = offset + 2
        dst = offset + 9
        for col_idx, col_letter in enumerate(("A", "B", "C", "D", "E", "F", "G"), start=10):
            ws.cell(dst, col_idx, f"='PK_Sample_Delta'!{col_letter}{src}")
    _style_table(ws, 8, 8 + max(min(assay_count, 12), 1), 10, 16)

    ws["A15"] = "Assay Watchlist"
    ws["A15"].font = Font(size=13, bold=True)
    watch_headers = ["Assay", "Files", "PatientFiles", "ControlFiles", "LadderReview", "AvgR2", "PartialFits", "ReviewRate"]
    for idx, header in enumerate(watch_headers, start=1):
        ws.cell(16, idx, header)
    for offset in range(0, min(assay_count, 12)):
        src = offset + 2
        dst = offset + 16
        for col_idx, col_letter in enumerate(("A", "B", "C", "D", "E", "F", "G", "H"), start=1):
            ws.cell(dst, col_idx, f"='Assay_Summary'!{col_letter}{src}")
    _style_table(ws, 16, 16 + max(min(assay_count, 12), 1), 1, 8)

    ws.conditional_formatting.add(
        f"E17:E{16 + max(min(assay_count, 12), 1)}",
        ColorScaleRule(start_type="num", start_value=0, start_color="E2F0D9", mid_type="percentile", mid_value=50, mid_color="FFE699", end_type="max", end_color="F4CCCC"),
    )
    ws.freeze_panes = "A8"


def _add_dashboard_charts(ws, assay_ws, pk_sample_ws) -> None:
    status_chart = BarChart()
    status_chart.title = "Ladder QC Status"
    status_chart.y_axis.title = "Runs"
    status_chart.height = 7
    status_chart.width = 8
    status_chart.add_data(Reference(ws, min_col=3, min_row=8, max_row=11), titles_from_data=True)
    status_chart.set_categories(Reference(ws, min_col=1, min_row=9, max_row=11))
    ws.add_chart(status_chart, "A30")

    assay_chart = BarChart()
    assay_chart.title = "Files by Assay"
    assay_chart.y_axis.title = "Files"
    assay_chart.height = 7
    assay_chart.width = 8
    max_row = max(assay_ws.max_row, 2)
    assay_chart.add_data(Reference(assay_ws, min_col=2, min_row=1, max_row=max_row), titles_from_data=True)
    assay_chart.set_categories(Reference(assay_ws, min_col=1, min_row=2, max_row=max_row))
    ws.add_chart(assay_chart, "E30")

    pk_chart = BarChart()
    pk_chart.title = "PK Sample Mean |delta bp| by Assay"
    pk_chart.y_axis.title = "|delta bp|"
    pk_chart.height = 7
    pk_chart.width = 8
    max_row = max(pk_sample_ws.max_row, 2)
    pk_chart.add_data(Reference(pk_sample_ws, min_col=3, min_row=1, max_row=max_row), titles_from_data=True)
    pk_chart.set_categories(Reference(pk_sample_ws, min_col=1, min_row=2, max_row=max_row))
    ws.add_chart(pk_chart, "J30")

