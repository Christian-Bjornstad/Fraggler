# Fraggler Diagnostics v1.2.0

## Highlights
- Standardized Fraggler as a native **PyQt desktop app** across macOS, Linux, and Windows.
- Reworked packaging around a single shared build contract so all desktop artifacts are produced from the same pipeline.
- Added a clearer **offline-first Linux delivery path** for the Fedora 35 work machine.

## What Changed

### Unified Desktop Packaging
- `qt_app.py` is now the canonical packaged application entrypoint.
- Build scripts for macOS, Linux, and Windows now route through the same desktop packaging logic.
- Release outputs are normalized so GitHub assets can follow a predictable naming scheme.

### Linux Offline Bundle
- Linux packaging is now designed around a **portable offline Fedora 35 bundle**.
- The Linux artifact is intended for:
  - `x86_64`
  - `glibc >= 2.31`
  - `QT_QPA_PLATFORM=xcb`
- Critical Qt/XCB-era runtime libraries are bundled into the offline release instead of assuming they already exist on the host.
- Deployment docs were updated so the Linux work machine setup is documented directly in the main branch.

### Runtime Behavior
- Packaged desktop builds now default to native desktop startup behavior.
- Legacy embedded/local Panel server behavior is no longer the default packaged mode.
- macOS keeps the `qt.conf` bundle fix to reduce translocation-related launch issues.

### Packaging Docs and Release Layout
- Packaging docs now describe the real desktop product instead of the older browser-first packaging model.
- Expected release assets are now:
  - `Fraggler_macOS.zip`
  - `Fraggler_Linux_offline.zip`
  - `Fraggler_Windows.zip`

## Recommended GitHub Release Assets
- `Fraggler_macOS.zip`
- `Fraggler_Linux_offline.zip`
- `Fraggler_Windows.zip`

## Validation Performed
- Unit test suite passes from the repo:
  - `python3 -m unittest discover -s tests -q`
- Verified the updated macOS packaging flow and artifact layout.
- Confirmed the new macOS release zip is produced under `dist/releases/`.

## Notes
- Linux runtime validation on the actual Fedora 35 machine is still the final real-world check for the offline bundle.
- macOS builds remain unsigned unless separately signed/notarized before release.
- PyInstaller still warns about `libomp.dylib` from `sklearn` during macOS packaging; it did not block the build, but it is worth keeping in mind for hardened-runtime/signing work later.
