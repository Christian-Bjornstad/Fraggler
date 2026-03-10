#!/bin/bash
# System dependencies for building Fraggler on Linux (offline)
# Run this as: sudo apt-get install ... (or on air-gapped: copy packages manually)

echo "=========================================="
echo "  Linux System Dependencies for Fraggler"
echo "=========================================="

# Core build tools
echo "Installing build tools..."
apt-get install -y gcc g++ make wget curl

# Qt6 and X11 dependencies (required for PyQt6)
echo "Installing Qt6/X11 dependencies..."
apt-get install -y \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    libx11-xcb1 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-xinerama0 \
    libxcb-xfixes0 \
    libxkbcommon-x11-0 \
    libxcb-shape0 \
    libxcb-glx0 \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libdbus-1-3 \
    libfontconfig1 \
    libfreetype6 \
    libxcb1 \
    libxcomposite1 \
    libxdamage1 \
    libxrandr2 \
    libxtst6 \
    libpango-1.0-0 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libcups2 \
    libdrm2 \
    libgbm1 \
    libasound2

echo ""
echo "Done! Now install Python dependencies:"
echo "  pip install --no-index --find-links=linux_offline_deps/ -r requirements.txt"
echo ""
echo "Then build:"
echo "  python build_qt.py"
