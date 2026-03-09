"""
Fraggler Diagnostics — App Entry Point (templateless)

panel serve app.py --port 5078 --allow-websocket-origin=localhost:5078
"""
from pathlib import Path
import panel as pn

_CSS_PATH = Path(__file__).resolve().parent / "assets" / "app.css"
_APP_CSS = _CSS_PATH.read_text(encoding="utf-8") if _CSS_PATH.exists() else ""

pn.extension(
    "plotly",
    "tabulator",
    sizing_mode="stretch_width",
    notifications=True,
    loading_spinner="dots",
    loading_color="#2563eb",  # Blue
    raw_css=[_APP_CSS],
)

from gui.main import build_app

build_app().servable(title="Fraggler Diagnostics")
