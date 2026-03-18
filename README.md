# Fraggler Diagnostics
A modern, standalone clinical diagnostics tool for automated processing, visualization, and quality control (QC) of DNA fragment analysis data (FSA traces). Developed with PyQt6 and optimized for capillary electrophoresis pipelines.

![App Icon](assets/app_icon.png)

## Features
- **Cross-Platform**: Run natively on macOS, Windows, or Linux via compiled binary executables.
- **Automated Processing**: Point the application at your FSA data folders; it automatically detects Patient sample sets and QC runs, processes the peaks, aligns scaling sizes, and determines pass/fail quality ratings natively.
- **Unified Batch Execution**: Add multiple source directories at once. Fraggler automatically builds nested pipeline jobs and isolates tracking parameters.
- **Clinical DIT Reports**: Generates interactive HTML graphical documents and aggregates QC logs into trending Excel spread-sheets to map ladder drift over time.
- **Secure/Offline**: No data leaves the application; all Plotly/HTML reporting tools and graphical assets are bundled natively within the app.

---

## Download & Installation

You can download the latest pre-compiled executables from the **[Releases](https://github.com/Christian-Bjornstad/Fraggler/releases)** page on GitHub! 

1. **macOS**: Download the `Fraggler.app` bundle and double-click to run. (You may need to bypass Apple's unsigned app warning by right-clicking it and selecting `Open`).
2. **Windows**: Download `Fraggler.exe` and execute it directly.
3. **Linux**: Download the ELF executable and run it via terminal (`./Fraggler`).

---

## Building from Source

If you wish to compile the standalone executables yourself, clone this repository and use the provided automated build scripts! 

#### Prerequisites:
- Python 3.10+
- macOS (For native `.app` builds)
- Docker Desktop (For cross-compiling the Windows and Linux binaries from your Mac)
- Build dependencies are pinned in `packaging/build-requirements.txt`

#### 1. macOS Build
Runs PyInstaller directly against your local Python environment to generate a `.app` bundle.
```bash
./packaging/build_mac.sh
```
*Wait for completion. You will find `Fraggler.app` located in the `/dist` directory.*

#### 2. Windows Build (via Docker)
Spins up an isolated Wine container (`tobix/pywine:3.10`) inside Docker to emulate a Windows system, installing PyInstaller and cross-compiling the Python codebase into a native Windows `.exe`.
```bash
./packaging/build_windows.sh
```
*Wait for completion. You will find `Fraggler.exe` exported to your local `/dist` directory.*

#### 3. Linux Build (via Docker)
Spins up an isolated Debian-slim Python Docker container, compiling the PyInstaller bundle into an ELF executable compatible with Ubuntu/Linux distributions.
```bash
./packaging/build_linux.sh
```
*Wait for completion. You will find the `Fraggler` linux executable exported to your local `/dist` directory.*

## Architecture
- **Language**: Python 3.10
- **GUI**: PyQt6 (Native integration, system tray, graphical routing)
- **Data Engineering**: Pandas, Numpy, Scikit-Learn (peak parsing and sequence alignments)
- **Visualization**: Plotly, Jinja2, HTML (interactive browser-facing diagnostics matrices)

## Troubleshooting
**"IndexError: tuple index out of range" during Build**
> During `PyInstaller` parsing, the build uses an internal monkey-patch embedded into `build_qt.py` to bypass a known Python 3.10.0 bytecode compilation `match` bug inside the overarching Plotly/Pandas dependencies. This acts as a silent failsafe that swallows the bug and permits the EXE to wrap correctly. 
