"""
Fraggler Diagnostics — Reports Tab (NEW)

Browse generated HTML reports, open inline, and print to PDF.
"""
from __future__ import annotations

import panel as pn
from pathlib import Path

from config import APP_SETTINGS
from gui.components import VSpace, section_header, make_card


def make_reports_tab() -> pn.Column:
    # ─── Controls ─────────────────────────────────────────────────────
    folder_input = pn.widgets.TextInput(
        name="Reports Folder",
        value=APP_SETTINGS.get("batch", {}).get("output_base", ""),
        sizing_mode="stretch_width",
        placeholder="/path/to/output/folder"
    )
    scan_btn = pn.widgets.Button(
        name="🔍  Scan for Reports", button_type="default", width=180, height=42,
        styles={"font-weight": "600"}
    )
    filter_input = pn.widgets.TextInput(
        placeholder="Filter by Patient ID / filename...",
        width=260,
        styles={"font-size": "13px"}
    )

    status_md = pn.pane.HTML(
        '<div style="color:#94a3b8; font-size:13px">Enter a folder path and click Scan.</div>',
        sizing_mode="stretch_width"
    )

    # File list
    file_select = pn.widgets.Select(
        name="Select Report",
        options=[],
        sizing_mode="stretch_width",
    )
    open_btn = pn.widgets.Button(
        name="📖 Open Report", button_type="primary", width=160, height=42,
        disabled=True, styles={"font-weight": "600"}
    )
    print_btn = pn.widgets.Button(
        name="🖨 Print / Save PDF", button_type="default", width=180, height=42,
        disabled=True, styles={"font-weight": "600"}
    )
    open_browser_btn = pn.widgets.Button(
        name="🌐 Open in Browser", button_type="default", width=180, height=42,
        disabled=True, styles={"font-weight": "600"}
    )

    # Viewer
    viewer = pn.pane.HTML(
        '<div style="display:flex; align-items:center; justify-content:center; height:500px; color:#475569; font-size:14px">No report loaded. Scan a folder and select a report.</div>',
        sizing_mode="stretch_both",
        min_height=600,
        styles={"border": "1px solid #2e3650", "border-radius": "8px", "overflow": "hidden"}
    )

    _all_paths: dict[str, Path] = {}
    _filtered_paths: dict[str, Path] = {}

    # ─── Scan ─────────────────────────────────────────────────────────
    def on_scan(event):
        nonlocal _all_paths, _filtered_paths
        folder = Path(folder_input.value).expanduser() if folder_input.value else None
        if not folder or not folder.exists():
            status_md.object = '<div style="color:#ef4444; font-size:13px">❌ Invalid folder path.</div>'
            return

        html_files = sorted(folder.rglob("*.html"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not html_files:
            status_md.object = '<div style="color:#f59e0b; font-size:13px">⚠️ No HTML reports found in this folder.</div>'
            _all_paths = {}
            _filtered_paths = {}
            file_select.options = []
            open_btn.disabled = True
            print_btn.disabled = True
            open_browser_btn.disabled = True
            return

        _all_paths = {f"{p.parent.name} / {p.name}": p for p in html_files}
        status_md.object = f'<div style="color:#22c55e; font-size:13px">✅ Found <strong>{len(html_files)}</strong> HTML reports.</div>'
        _apply_filter(filter_input.value)

    scan_btn.on_click(on_scan)

    # ─── Filter ───────────────────────────────────────────────────────
    def _apply_filter(search: str = ""):
        nonlocal _filtered_paths
        if search:
            _filtered_paths = {k: v for k, v in _all_paths.items() if search.lower() in k.lower()}
        else:
            _filtered_paths = dict(_all_paths)
        file_select.options = list(_filtered_paths.keys()) if _filtered_paths else []
        open_btn.disabled = not bool(_filtered_paths)
        print_btn.disabled = not bool(_filtered_paths)
        open_browser_btn.disabled = not bool(_filtered_paths)

    def _on_filter(event):
        _apply_filter(event.new)

    filter_input.param.watch(_on_filter, "value")

    # ─── Open ─────────────────────────────────────────────────────────
    def on_open(event):
        key = file_select.value
        if not key or key not in _filtered_paths:
            return
        path = _filtered_paths[key]
        try:
            content = path.read_text(encoding="utf-8", errors="replace")
            # Embed the report directly in an iframe
            escaped = content.replace("'", "&#39;")
            viewer.object = f"""<iframe
  srcdoc='{escaped}'
  id="report-iframe"
  style="width:100%; height:700px; border:none; border-radius:8px"
></iframe>"""
            print_btn.disabled = False
        except Exception as e:
            viewer.object = f'<div style="color:#ef4444; padding:20px">Error loading report: {e}</div>'

    open_btn.on_click(on_open)

    # Auto-open when selection changes
    def _on_select(event):
        if event.new and _filtered_paths:
            on_open(None)

    file_select.param.watch(_on_select, "value")

    # ─── Print ────────────────────────────────────────────────────────
    _print_js = pn.pane.HTML(
        '<script id="print-trigger"></script>',
        width=0, height=0
    )

    def on_print(event):
        # Trigger print on the embedded iframe content
        _print_js.object = """<script>
(function(){
  var iframe = document.getElementById('report-iframe');
  if(iframe) {
    iframe.contentWindow.focus();
    iframe.contentWindow.print();
  } else {
    window.print();
  }
})();
</script>"""

    print_btn.on_click(on_print)

    # ─── Open in browser ──────────────────────────────────────────────
    def on_open_browser(event):
        import subprocess, sys
        key = file_select.value
        if not key or key not in _filtered_paths:
            return
        path = _filtered_paths[key]
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])
        elif sys.platform == "win32":
            subprocess.Popen(["explorer", str(path)])
        else:
            subprocess.Popen(["xdg-open", str(path)])

    open_browser_btn.on_click(on_open_browser)

    # ─── Layout ───────────────────────────────────────────────────────
    return pn.Column(
        section_header("Reports Viewer", "Browse, preview and print HTML reports from batch runs"),
        VSpace(8),

        pn.Row(
            folder_input,
            scan_btn,
            styles={"gap": "10px", "align-items": "end"}
        ),
        status_md,
        VSpace(10),

        pn.Row(
            file_select,
            filter_input,
            styles={"gap": "12px", "align-items": "end"}
        ),
        VSpace(6),

        pn.Row(
            open_btn, open_browser_btn, print_btn,
            _print_js,
            styles={"gap": "10px", "align-items": "center", "flex-wrap": "wrap"}
        ),
        VSpace(10),

        pn.pane.HTML(
            '<div style="font-size:11px; font-weight:600; color:#94a3b8; text-transform:uppercase; letter-spacing:0.8px; margin-bottom:8px">Report Preview</div>',
            sizing_mode="stretch_width"
        ),
        viewer,

        sizing_mode="stretch_both",
        styles={"padding": "20px 24px", "gap": "0", "max-width": "1300px"},
    )
