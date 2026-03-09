# PyInstaller hook for Panel
# Ensures Panel template files and static resources are included.

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = collect_data_files("panel")
hiddenimports = collect_submodules("panel")
