# PyInstaller hook for fraggler
# Ensures the fraggler package (ladders data) is included.

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = collect_data_files("fraggler")
hiddenimports = collect_submodules("fraggler")
