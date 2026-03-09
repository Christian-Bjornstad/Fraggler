# Fraggler Diagnostics — Packaging

Build standalone executables for **macOS**, **Linux**, and **Windows**.

## Quick Start

### Build for Mac (current machine)
```bash
cd /Users/christian/Desktop/OUS
chmod +x packaging/build_mac.sh
./packaging/build_mac.sh
```

Result: `packaging/dist/fraggler-diagnostics/fraggler-diagnostics`

### Build for Linux (via Docker)
```bash
cd /Users/christian/Desktop/OUS
chmod +x packaging/build_linux.sh
./packaging/build_linux.sh
```

Result: `packaging/dist/fraggler-diagnostics-linux/fraggler-diagnostics`

### Build for Windows
Run on a Windows machine with Python 3.10+:
```cmd
packaging\build_windows.bat
```

## Running the Executable

```bash
# Mac / Linux
./fraggler-diagnostics

# Windows
fraggler-diagnostics.exe
```

The app will:
1. Start a local web server on port 5078
2. Automatically open your browser at `http://localhost:5078/app`
3. Show logs in the terminal

Press `Ctrl+C` in the terminal to stop.

## How It Works

```
┌─────────────────────────────────────────────┐
│  fraggler-diagnostics (executable)          │
│  ┌───────────────────┐                      │
│  │  launcher.py       │ ← Entry point       │
│  │  - Starts Panel    │                      │
│  │  - Opens browser   │                      │
│  └────────┬──────────┘                      │
│           ↓                                  │
│  ┌───────────────────┐                      │
│  │  _internal/        │ ← Bundled deps       │
│  │  - Python 3.10     │                      │
│  │  - numpy, scipy    │                      │
│  │  - panel, bokeh    │                      │
│  │  - plotly.js       │                      │
│  │  - all app code    │                      │
│  └───────────────────┘                      │
└─────────────────────────────────────────────┘
```

## Files

| File | Purpose |
|------|---------|
| `launcher.py` | Entry point — starts server & opens browser |
| `fraggler_diagnostics.spec` | PyInstaller config — what to bundle |
| `hooks/hook-*.py` | PyInstaller hooks for Panel/Bokeh/fraggler |
| `build_mac.sh` | Mac build script |
| `build_linux.sh` | Linux build via Docker |
| `build_windows.bat` | Windows build script |
| `Dockerfile.linux` | Docker image for Linux builds |

## Troubleshooting

### "Port 5078 already in use"
Another instance is running. Kill it:
```bash
# Mac/Linux
lsof -i :5078 -t | xargs kill

# Windows
netstat -ano | findstr :5078
taskkill /PID <pid> /F
```

### Large executable size (~150-250 MB)
This is expected — it includes Python + numpy + scipy + sklearn + pandas + plotly.js.

### Linux executable won't run
Make sure it's executable:
```bash
chmod +x fraggler-diagnostics
```
