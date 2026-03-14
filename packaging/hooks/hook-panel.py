# PyInstaller hook for Panel
# Ensures Panel template files and static resources are included.

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

datas = collect_data_files("panel")
for package_name in ("panel", "param", "pyviz_comms"):
    datas += copy_metadata(package_name)
hiddenimports = collect_submodules("panel")
