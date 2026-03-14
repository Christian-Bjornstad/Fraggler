# Fraggler Diagnostics — Linux Offline Guide

This guide is for the **Fedora 35 x86_64 work machine** and other Linux systems with similar compatibility characteristics.

Fraggler is shipped on Linux as a **portable offline desktop bundle**.

## Recommended Artifact

Use:
- `Fraggler_Linux_offline.zip`

The bundle includes:
- the `Fraggler` launcher
- the `_internal` runtime folder from PyInstaller
- bundled Qt/XCB runtime libraries needed for the Fedora 35 target
- a Linux-specific README inside the bundle

## Supported Runtime Assumptions

- `glibc >= 2.31`
- `x86_64`
- offline deployment supported
- packaged Linux runs with `QT_QPA_PLATFORM=xcb`

Fedora 35 typically reports glibc 2.34, so it satisfies the baseline.

## Install and Run

```bash
unzip Fraggler_Linux_offline.zip -d ~/Fraggler
cd ~/Fraggler/Fraggler_Linux
chmod +x Fraggler
./Fraggler
```

## Validation Checklist

Check glibc:
```bash
ldd --version
```

Check runtime dependencies:
```bash
ldd ./Fraggler
```

Expected:
- no unresolved XCB/Qt runtime libraries that were meant to be bundled
- the application stays in X11/xcb mode

## Troubleshooting

### "Symbol lookup error" or "Library not found"
- Confirm the full zip was extracted
- Keep the `_internal` folder next to `Fraggler`
- Run `ldd ./Fraggler` to identify unresolved libraries

### Wayland vs X11
Fraggler forces **X11/xcb** in packaged Linux builds for compatibility with the Fedora 35 target environment.

### Permission Denied
Run:
```bash
chmod +x Fraggler
```
