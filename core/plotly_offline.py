"""
Offline Plotly script tag provider.

Reads plotly-3.1.0.min.js from the bundled assets/ directory and returns
a <script> tag with the FULL JS inlined.  No internet required.  No external
file needs to be shipped alongside the HTML.
"""
from __future__ import annotations

from pathlib import Path

_ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"
_PLOTLY_JS = _ASSETS_DIR / "plotly-3.1.0.min.js"

# Cache the JS content so we only read it once per process
_PLOTLY_INLINE_CACHE: str | None = None


def plotly_inline_script_tag() -> str:
    """
    Returns a ``<script>…</script>`` tag with the full plotly.js inlined.

    The JS is read from ``assets/plotly-3.1.0.min.js`` on first call and
    cached for subsequent calls.
    """
    global _PLOTLY_INLINE_CACHE
    if _PLOTLY_INLINE_CACHE is None:
        if not _PLOTLY_JS.exists():
            raise FileNotFoundError(
                f"Bundled plotly JS not found at {_PLOTLY_JS}.  "
                f"Please ensure assets/plotly-3.1.0.min.js exists."
            )
        # Use errors="replace" to avoid crashing on unexpected characters in minified JS
        _PLOTLY_INLINE_CACHE = _PLOTLY_JS.read_text(encoding="utf-8", errors="replace")
    return f"<script>{_PLOTLY_INLINE_CACHE}</script>"


# ------------------------------------------------------------------
# Backward-compatible aliases
# ------------------------------------------------------------------
def local_plotly_tag(out_dir: "Path | None" = None, version: str = "3.1.0") -> str:
    """Drop-in replacement for the old ``plotly_local.local_plotly_tag``."""
    return plotly_inline_script_tag()
