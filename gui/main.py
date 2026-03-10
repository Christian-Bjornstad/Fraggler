"""
Fraggler Diagnostics — Main Layout (templateless, light software theme)

Modern software look: 
- Dark navy sidebar
- Light clean background
- White cards with subtle borders
- Blue accent colors
- Hides default Panel header/theme toggle via CSS
"""
from __future__ import annotations
from pathlib import Path
import panel as pn

_CSS_PATH = Path(__file__).resolve().parent.parent / "assets" / "app.css"
_APP_CSS = _CSS_PATH.read_text(encoding="utf-8") if _CSS_PATH.exists() else ""

def build_app() -> pn.Column:
    """Build the full application as a pn.Column (templateless)."""
    from gui.tab_batch    import make_batch_tab
    from gui.tab_qc       import make_qc_tab
    from gui.tab_log      import make_log_tab
    from gui.tab_settings import make_settings_tab
    from gui.analysis_selector import make_analysis_selector

    def _wrap(content: pn.viewable.Viewable) -> pn.Column:
        return pn.Column(
            content, 
            sizing_mode="stretch_both",
            styles={"background": "transparent"}
        )

    tabs = pn.Tabs(
        ("Batch",    _wrap(make_batch_tab())),
        ("QC",       _wrap(make_qc_tab())),
        ("Log",      _wrap(make_log_tab())),
        ("Settings", _wrap(make_settings_tab())),
        tabs_location="left",
        sizing_mode="stretch_both",
        dynamic=False,
        css_classes=["main-tabs"],
        stylesheets=[_APP_CSS],
        styles={"background": "transparent"},
    )

    header = pn.Row(
        pn.pane.HTML('<div style="font-size:20px; font-weight:bold; color:var(--accent)">Fraggler Diagnostics</div>'),
        pn.HSpacer(),
        pn.Column(make_analysis_selector(), width=250),
        sizing_mode="stretch_width",
        styles={"padding": "10px 20px", "background": "var(--bg-card)", "border-bottom": "1px solid var(--border)"},
        css_classes=["app-header"]
    )

    return pn.Column(
        header,
        tabs,
        sizing_mode="stretch_both",
        styles={"background": "var(--bg-app)"}
    )
