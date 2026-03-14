# PyInstaller hook for Bokeh
# Ensures Bokeh model definitions and static resources are included.

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

datas = collect_data_files("bokeh")
datas += copy_metadata("bokeh")
hiddenimports = collect_submodules("bokeh")
