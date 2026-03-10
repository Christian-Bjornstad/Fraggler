"""
Analysis Selector Component.
"""
from __future__ import annotations
import panel as pn
from config import APP_SETTINGS, save_settings

def make_analysis_selector() -> pn.widgets.Select:
    active = APP_SETTINGS.get("active_analysis", "clonality")
    
    selector = pn.widgets.Select(
        name="Select Analysis",
        options={
            "Klonalitet (Clonality)": "clonality",
            "FLT3 (Skeleton)": "flt3"
        },
        value=active,
        sizing_mode="stretch_width",
        css_classes=["analysis-selector"]
    )
    
    def on_change(event):
        APP_SETTINGS["active_analysis"] = event.new
        save_settings(APP_SETTINGS)
        # We might need to refresh the page or update components
        pn.state.notifications.info(f"Analysis switched to: {event.new.upper()}", duration=3000)
        # Note: True "refresh" of components might require a callback 
        # to the parent to rebuild or reload.
        # For now, we rely on the fact that registry.py pulls from APP_SETTINGS dynamically.
        
    selector.param.watch(on_change, "value")
    return selector
