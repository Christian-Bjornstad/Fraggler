# Fraggler Diagnostics — Linux (Offline Guide)

This guide is for running Fraggler on a Linux machine (like Fedora 35) without an internet connection.

## 🚀 Recommended: Full Portable Bundle

We have created a "Full Bundle" that includes its own copy of critical system libraries (like `libxcb-cursor`) so you don't need to install anything on the host.

1. **Locate the package**: Use `Fraggler_Linux_offline.zip`.
2. **Transfer to Linux**: Use a USB drive to copy it to your Fedora machine.
3. **Extract and Run**:
   ```bash
   unzip Fraggler_Linux_offline.zip -d ~/Fraggler
   cd ~/Fraggler/dist/Fraggler_Linux
   chmod +x Fraggler
   ./Fraggler
   ```

## 🛠 Prerequisites

- **GLIBC 2.31 or newer**: Check with `ldd --version` (Fedora 35 is 2.34, so it works).
- **Architecture**: 64-bit Intel/AMD (x86_64).

## ❓ Troubleshooting

### "Symbol lookup error" or "Library not found"
The **v2 Offline** version includes all known missing libraries. If you still see errors, check that you extracted the entire zip file, as the `_internal` folder contains the required `.so` files.

### Wayland vs X11
The application is configured to force the **X11 (xcb)** platform for maximum compatibility. This is handled automatically by the `Fraggler` binary.

### Permission Denied
Run `chmod +x Fraggler` in the terminal inside the `Fraggler_Linux` folder.
