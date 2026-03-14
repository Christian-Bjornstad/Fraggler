from __future__ import annotations

import dis
import os
import shutil
import sys
import zipfile
from pathlib import Path

import PyInstaller.__main__

# Monkey-patch dis._get_const_info to swallow known Python 3.10 bytecode IndexErrors
_orig_get_const_info = getattr(dis, '_get_const_info', None)
if _orig_get_const_info:
    def _patched_get_const_info(arg, constants):
        try:
            return _orig_get_const_info(arg, constants)
        except IndexError:
            return arg, repr(arg)
    dis._get_const_info = _patched_get_const_info

APP_NAME = "Fraggler"
BUNDLE_ID = "com.christian-bjornstad.fraggler"
PROJECT_ROOT = Path(__file__).resolve().parent
DIST_DIR = PROJECT_ROOT / "dist"
RELEASE_DIR = DIST_DIR / "releases"
HOOK_DIR = PROJECT_ROOT / "packaging" / "hooks"
LEGACY_LINUX_GUIDE = PROJECT_ROOT / "LINUX_GUIDE.md"

COMMON_HIDDEN_IMPORTS = [
    "PyQt6",
    "pandas",
    "plotly",
    "core.analyses.clonality.config",
    "core.analyses.clonality.classification",
    "core.analyses.clonality.pipeline",
    "core.analyses.flt3.config",
    "core.analyses.flt3.classification",
    "core.analyses.flt3.pipeline",
]

COMMON_DATAS = [
    ("assets", "assets"),
]

LINUX_BINARY_CANDIDATES = [
    "/usr/lib/x86_64-linux-gnu/libxcb-cursor.so.0",
    "/usr/lib/x86_64-linux-gnu/libxcb-icccm.so.4",
    "/usr/lib/x86_64-linux-gnu/libxcb-image.so.0",
    "/usr/lib/x86_64-linux-gnu/libxcb-keysyms.so.1",
    "/usr/lib/x86_64-linux-gnu/libxcb-randr.so.0",
    "/usr/lib/x86_64-linux-gnu/libxcb-render-util.so.0",
    "/usr/lib/x86_64-linux-gnu/libxcb-xinerama.so.0",
    "/usr/lib/x86_64-linux-gnu/libxcb-xfixes.so.0",
    "/usr/lib/x86_64-linux-gnu/libxkbcommon-x11.so.0",
    "/usr/lib/x86_64-linux-gnu/libdbus-1.so.3",
    "/usr/lib/x86_64-linux-gnu/libfontconfig.so.1",
    "/usr/lib/x86_64-linux-gnu/libGL.so.1",
    "/usr/lib/x86_64-linux-gnu/libEGL.so.1",
    "/usr/lib/x86_64-linux-gnu/libglib-2.0.so.0",
]


def _pyinstaller_sep() -> str:
    return ";" if sys.platform == "win32" else ":"


def _format_data_arg(src: str, dest: str) -> str:
    return f"--add-data={src}{_pyinstaller_sep()}{dest}"


def _format_binary_arg(src: str, dest: str) -> str:
    return f"--add-binary={src}{_pyinstaller_sep()}{dest}"


def _collect_linux_binaries() -> list[str]:
    args: list[str] = []
    for path_str in LINUX_BINARY_CANDIDATES:
        path = Path(path_str)
        if path.exists():
            print(f"Bundling Linux runtime library: {path}")
            args.append(_format_binary_arg(str(path), "."))
    return args


def _build_pyinstaller_args() -> list[str]:
    args = [
        "qt_app.py",
        f"--name={APP_NAME}",
        "--noconfirm",
        "--clean",
        "--windowed",
        f"--distpath={DIST_DIR}",
        f"--workpath={PROJECT_ROOT / 'build'}",
        f"--specpath={PROJECT_ROOT}",
        f"--additional-hooks-dir={HOOK_DIR}",
        f"--runtime-hook={HOOK_DIR / 'runtime_desktop.py'}",
    ]

    for src, dest in COMMON_DATAS:
        args.append(_format_data_arg(src, dest))
    for mod in COMMON_HIDDEN_IMPORTS:
        args.append(f"--hidden-import={mod}")

    if sys.platform == "darwin":
        args.append("--icon=assets/app_icon.icns")
        args.append(f"--osx-bundle-identifier={BUNDLE_ID}")
    elif sys.platform == "win32":
        args.append("--icon=assets/app_icon.ico")
    elif sys.platform == "linux":
        args.extend(_collect_linux_binaries())

    return args


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _zip_path(src: Path, zip_path: Path, root_name: str | None = None) -> Path:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        if src.is_dir():
            for file_path in sorted(src.rglob("*")):
                if file_path.is_dir():
                    continue
                rel = file_path.relative_to(src)
                arcname = Path(root_name or src.name) / rel
                zf.write(file_path, arcname.as_posix())
        else:
            zf.write(src, (root_name or src.name))
    return zip_path


def _post_build_mac() -> None:
    resources_dir = DIST_DIR / f"{APP_NAME}.app" / "Contents" / "Resources"
    if resources_dir.exists():
        qt_conf_path = resources_dir / "qt.conf"
        print(f"Creating {qt_conf_path} to fix translocation crashes...")
        _write_text(qt_conf_path, "[Paths]\nPrefix = .\n")


def _linux_readme_text() -> str:
    return """Fraggler Diagnostics — Fedora 35 Offline Bundle

Target:
- Fedora 35 (or other x86_64 Linux with glibc 2.31+)
- Offline deployment supported
- Native PyQt desktop application

Run:
1. Extract this zip fully.
2. Open a terminal in the extracted Fraggler_Linux folder.
3. Make the launcher executable:
   chmod +x Fraggler
4. Start the app:
   ./Fraggler

Compatibility:
- The app forces QT_QPA_PLATFORM=xcb for X11 compatibility.
- Critical XCB/Qt runtime libraries are bundled in this portable folder.

Validation checklist:
- ldd --version reports glibc 2.31 or newer.
- ldd ./Fraggler does not show missing bundled runtime dependencies.
- The whole extracted folder is kept intact, including the _internal directory.

Troubleshooting:
- If launch fails, run:
  ldd ./Fraggler
- If you see missing library errors, confirm the full zip was extracted.
"""


def _stage_release_dir(name: str) -> Path:
    target = DIST_DIR / name
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    return target


def _post_build_linux() -> None:
    source_dir = DIST_DIR / APP_NAME
    release_dir = _stage_release_dir("Fraggler_Linux")
    shutil.copytree(source_dir, release_dir, dirs_exist_ok=True)

    if LEGACY_LINUX_GUIDE.exists():
        shutil.copy2(LEGACY_LINUX_GUIDE, release_dir / "LINUX_GUIDE.md")
    _write_text(release_dir / "README.txt", _linux_readme_text())

    launcher = release_dir / APP_NAME
    if launcher.exists():
        launcher.chmod(0o755)

    zip_path = RELEASE_DIR / "Fraggler_Linux_offline.zip"
    _zip_path(release_dir, zip_path, root_name=release_dir.name)
    print(f"Created Linux offline bundle: {zip_path}")


def _post_build_windows() -> None:
    source_dir = DIST_DIR / APP_NAME
    release_dir = _stage_release_dir("Fraggler_Windows")
    shutil.copytree(source_dir, release_dir, dirs_exist_ok=True)
    zip_path = RELEASE_DIR / "Fraggler_Windows.zip"
    _zip_path(release_dir, zip_path, root_name=release_dir.name)
    print(f"Created Windows bundle: {zip_path}")


def _post_build_generic_desktop() -> None:
    app_path = DIST_DIR / f"{APP_NAME}.app"
    zip_path = RELEASE_DIR / "Fraggler_macOS.zip"
    _zip_path(app_path, zip_path, root_name=app_path.name)
    print(f"Created macOS bundle: {zip_path}")


def build_app() -> None:
    print(f"Building Fraggler Diagnostics desktop app for {sys.platform}...")
    PyInstaller.__main__.run(_build_pyinstaller_args())

    RELEASE_DIR.mkdir(parents=True, exist_ok=True)
    if sys.platform == "darwin":
        _post_build_mac()
        _post_build_generic_desktop()
    elif sys.platform == "linux":
        _post_build_linux()
    elif sys.platform == "win32":
        _post_build_windows()

    print("\nBuild complete! Check the 'dist' directory.")

if __name__ == "__main__":
    build_app()
