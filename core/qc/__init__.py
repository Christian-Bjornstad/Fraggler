"""
HemaFrag QC subpackage re-exports.
"""
from core.qc.qc_rules import QCRules, ASSAY_ALIASES_QC, normalize_assay_qc  # noqa


def __getattr__(name: str):
    if name == "build_qc_html":
        from core.qc.qc_html import build_qc_html
        return build_qc_html
    if name in {"update_excel_trends", "apply_pk_excel_styling"}:
        from core.qc.qc_excel import update_excel_trends, apply_pk_excel_styling
        return {"update_excel_trends": update_excel_trends, "apply_pk_excel_styling": apply_pk_excel_styling}[name]
    if name == "main":
        from core.qc.qc_main import main
        return main
    raise AttributeError(f"module {__name__} has no attribute {name}")
