"""
Fraggler Diagnostics — GUI Shared Components

Reusable widgets, badges, and layout helpers.
"""
from __future__ import annotations

import panel as pn

# Layout constants
GAP = "12px"

def VSpace(height: int = 10) -> pn.layout.Spacer:
    return pn.layout.Spacer(height=height, sizing_mode="stretch_width")

def HSpace(width: int = 10) -> pn.layout.Spacer:
    return pn.layout.Spacer(width=width, height=1)

def make_card(title: str, *objects, collapsed: bool = False, **kwargs) -> pn.Card:
    styles = {"gap": GAP, "padding": "16px"}
    if "styles" in kwargs:
        styles.update(kwargs.pop("styles"))
    return pn.Card(
        *objects,
        title=title,
        collapsed=collapsed,
        sizing_mode="stretch_width",
        styles=styles,
        **kwargs
    )

def stat_card(value: str, label: str) -> pn.pane.HTML:
    """Small stat counter card."""
    return pn.pane.HTML(
        f"""<div class="stat-card">
  <div class="stat-value">{value}</div>
  <div class="stat-label">{label}</div>
</div>""",
        width=130,
    )

def badge(text: str, kind: str = "pending") -> str:
    """Return an HTML badge string. kind: pending|running|done|error|skipped"""
    return f'<span class="badge badge-{kind}">{text}</span>'

def section_header(title: str, subtitle: str = "") -> pn.pane.HTML:
    sub_html = f"<p style='color: var(--text-muted); font-size:13px; margin:4px 0 0'>{subtitle}</p>" if subtitle else ""
    return pn.pane.HTML(
        f"""<div style="margin-bottom: 4px;">
  <h1 style="font-size:22px; font-weight:700; color:var(--text); margin:0">{title}</h1>
  {sub_html}
</div>""",
        sizing_mode="stretch_width",
    )

def hero_row(*btns) -> pn.Row:
    """A row to hold the primary hero action buttons."""
    return pn.Row(
        *btns,
        styles={"gap": "12px", "align-items": "center", "flex-wrap": "wrap"},
    )
