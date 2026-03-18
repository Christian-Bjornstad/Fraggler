# Fraggler Diagnostics — Desktop Packaging

Fraggler is packaged as a **native PyQt desktop app** on all supported platforms.
The canonical packaged entrypoint is `qt_app.py`.

The embedded/local Panel server is now treated as legacy support only:
- it is not the default packaged startup path
- packaged builds disable it by default via a runtime hook
- it can still be re-enabled explicitly with `FRAGGLER_ENABLE_LEGACY_PANEL=1`

## Supported Release Outputs

Every platform is built from the same desktop packaging contract in `build_qt.py`.

Artifacts:
- macOS: `dist/Fraggler.app` and `dist/releases/Fraggler_macOS.zip`
- Linux: `dist/Fraggler_Linux` and `dist/releases/Fraggler_Linux_offline.zip`
- Windows: `dist/Fraggler_Windows` and `dist/releases/Fraggler_Windows.zip`

## Platform Targets

### macOS
- Native `.app` bundle
- Includes `qt.conf` post-build fix for macOS translocation compatibility

### Linux
- Primary target: **offline Fedora 35 x86_64 machine**
- Runtime assumption: **glibc 2.31 or newer**
- Desktop startup forces **`QT_QPA_PLATFORM=xcb`**
- Offline bundle includes critical XCB/Qt runtime libraries inside the artifact

### Windows
- Native desktop bundle with `Fraggler.exe`
- Same packaged app behavior as macOS/Linux

## Quick Start

### Build macOS
```bash
cd /path/to/OUS
./packaging/build_mac.sh
```

### Build Linux Offline Bundle
```bash
cd /path/to/OUS
./packaging/build_linux.sh
```

### Build Windows on Windows
```cmd
packaging\build_windows.bat
```

### Build Windows via Docker
```bash
cd /path/to/OUS
./packaging/build_windows.sh
```

Build-tool reproducibility:
- `packaging/build-requirements.txt` pins the packaging toolchain separately from runtime dependencies.
- The wrapper scripts and Docker build images install both `requirements.txt` and `packaging/build-requirements.txt`.
- Full CI/release automation is still out of scope; artifact launch validation remains a manual release check.

## Linux Offline Deployment

Use `dist/releases/Fraggler_Linux_offline.zip`.

Target machine:
- Fedora 35
- x86_64
- offline supported

Run on Linux:
```bash
unzip Fraggler_Linux_offline.zip -d ~/Fraggler
cd ~/Fraggler/Fraggler_Linux
chmod +x Fraggler
./Fraggler
```

Validation checklist:
```bash
ldd --version
ldd ./Fraggler
```

Expected:
- glibc 2.31+
- no unresolved bundled XCB/Qt runtime dependencies

## Files

| File | Purpose |
|------|---------|
| `build_qt.py` | Canonical PyInstaller desktop build contract |
| `packaging/build_mac.sh` | macOS wrapper around the shared desktop build |
| `packaging/build_linux.sh` | Docker-based Linux offline bundle build |
| `packaging/build_windows.bat` | Native Windows desktop build |
| `packaging/build_windows.sh` | Docker/Wine-based Windows desktop build |
| `packaging/hooks/runtime_desktop.py` | Packaged runtime defaults |
| `packaging/Dockerfile.linux` | Older-compatible Linux build image |
| `LINUX_GUIDE.md` | Fedora 35 deployment notes merged into mainline docs |

## Release Uploads

Recommended GitHub release assets:
- `Fraggler_macOS.zip`
- `Fraggler_Linux_offline.zip`
- `Fraggler_Windows.zip`

## Troubleshooting

### Linux launch fails with missing library errors
- Ensure the full zip was extracted, including `_internal`
- Run `ldd ./Fraggler`
- Keep the offline bundle contents together; do not move the binary out of the folder

### Linux opens with Wayland/X11 issues
Fraggler forces `QT_QPA_PLATFORM=xcb` in packaged Linux builds for compatibility.

### macOS unsigned app warning
Use right-click → `Open` the first time, or sign/notarize the zip before release.

### Windows or Linux build size is large
Expected. The bundle includes Python, scientific libraries, Qt, and packaged assets.
