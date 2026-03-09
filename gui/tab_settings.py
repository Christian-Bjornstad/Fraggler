"""
Fraggler Diagnostics — Settings Tab
"""
from __future__ import annotations

import panel as pn

from config import APP_SETTINGS, save_settings
from gui.components import make_card, VSpace


def make_settings_tab() -> pn.Column:
    status_md = pn.pane.HTML("", sizing_mode="stretch_width")

    # Global output default
    default_output = pn.widgets.TextInput(
        name="Global Default Output Folder",
        value=APP_SETTINGS.get("default_output", ""),
        sizing_mode="stretch_width",
        placeholder="/path/to/default/output"
    )
    save_output_btn = pn.widgets.Button(name="Save", button_type="default", width=100, height=38)

    def save_output(event):
        APP_SETTINGS["default_output"] = default_output.value
        save_settings(APP_SETTINGS)
        status_md.object = '<div style="color:var(--green); font-size:13px">✅ Default output folder saved.</div>'

    save_output_btn.on_click(save_output)

    # App info
    info_html = pn.pane.HTML("""
<div style="background:#ffffff; border:1px solid var(--border); border-radius:8px; padding:16px; font-size:13px; color:var(--text-dim); line-height:1.8">
  <div style="font-size:20px; font-weight:700; color:var(--text); margin-bottom:8px">Fraggler Diagnostics</div>
  <div><span style="color:var(--primary); font-weight:600">Version:</span> 2.0 — Professional Edition</div>
  <div><span style="color:var(--primary); font-weight:600">Framework:</span> Panel / Bokeh (localhost server)</div>
  <div><span style="color:var(--primary); font-weight:600">Analysis:</span> Fraggler library (fragment analysis)</div>
  <div style="margin-top:10px; font-size:11px; color:var(--muted)">Settings are stored in config.json alongside the application.</div>
</div>""", sizing_mode="stretch_width")

    return pn.Column(
        pn.pane.HTML('<div class="page-title">Settings</div><div class="page-sub">Application preferences and defaults</div>'),
        VSpace(8),

        make_card(
            "Default Paths",
            default_output,
            pn.Row(save_output_btn),
            status_md,
            css_classes=["fd-card"],
        ),
        VSpace(8),
        make_card("About", info_html, css_classes=["fd-card"]),

        sizing_mode="stretch_width",
        styles={"padding": "24px 28px", "gap": "0", "max-width": "800px"},
    )
