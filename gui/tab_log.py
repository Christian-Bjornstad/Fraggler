"""
Fraggler Diagnostics — Log Tab (Redesigned)

Terminal-style dark log with search filter and export.
"""
from __future__ import annotations

import panel as pn

from core.log import log_buffer
from gui.components import VSpace, section_header


def make_log_tab() -> pn.Column:
    auto_scroll = pn.widgets.Checkbox(name="Auto-scroll", value=True)

    def _render_log(text: str, search: str = "") -> str:
        lines = text.split("\n") if text else []
        if search:
            lines = [l for l in lines if search.lower() in l.lower()]
        html_lines = []
        for line in lines[-2000:]:  # cap at 2000 lines
            cls = ""
            if "[ERROR]" in line or "❌" in line:
                cls = "log-error"
            elif "[WARN]" in line or "⚠" in line:
                cls = "log-warn"
            elif "[INFO]" in line or "[BATCH]" in line or "✅" in line:
                cls = "log-info"
            elif "[SUCCESS]" in line:
                cls = "log-success"
            escaped = line.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            html_lines.append(f'<div class="{cls}">{escaped}</div>' if cls else f'<div>{escaped}</div>')

        scroll_js = """
<script>
(function(){
  var el = document.getElementById('log-terminal-inner');
  if(el) el.scrollTop = el.scrollHeight;
})();
</script>""" if auto_scroll.value else ""

        return f"""<div id="log-terminal-inner" class="log-terminal" style="min-height:480px; max-height:70vh; overflow-y:auto">
{''.join(html_lines) if html_lines else '<div style="color:#475569">No log entries yet.</div>'}
</div>{scroll_js}"""

    # Terminal-style HTML pane that mirrors log_buffer
    log_display = pn.pane.HTML(
        _render_log(""),
        sizing_mode="stretch_both",
        min_height=500,
    )

    search_input = pn.widgets.TextInput(
        placeholder="Filter log lines...",
        width=300,
        styles={"font-size": "13px"}
    )
    clear_btn = pn.widgets.Button(name="Clear", button_type="warning", width=100, height=36)
    export_btn = pn.widgets.Button(name="Export .txt", button_type="default", width=120, height=36)

    def _refresh(search: str = ""):
        log_display.object = _render_log(log_buffer.text, search)

    # Watch log buffer
    def _on_log_change(event):
        _refresh(search_input.value)

    log_buffer.param.watch(_on_log_change, "text")

    def _on_search(event):
        _refresh(event.new)

    search_input.param.watch(_on_search, "value")

    def on_clear(event):
        log_buffer.clear()
        log_display.object = _render_log("")

    def on_export(event):
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        path = f"/tmp/fraggler_log_{ts}.txt"
        with open(path, "w", encoding="utf-8") as f:
            f.write(log_buffer.text)
        log_buffer.text = log_buffer.text + f"\n[EXPORT] Log saved to {path}"

    clear_btn.on_click(on_clear)
    export_btn.on_click(on_export)

    return pn.Column(
        pn.pane.HTML('<div class="page-title">System Log</div><div class="page-sub">Real-time output from pipeline, QC, and batch jobs.</div>'),
        VSpace(8),
        pn.Row(
            search_input,
            pn.layout.Spacer(sizing_mode="stretch_width"),
            auto_scroll,
            clear_btn,
            export_btn,
            styles={"gap": "10px", "align-items": "center", "flex-wrap": "wrap"}
        ),
        VSpace(8),
        log_display,
        sizing_mode="stretch_both",
        styles={"padding": "24px 28px", "gap": "0", "max-width": "1300px"},
    )
