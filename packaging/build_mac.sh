#!/bin/bash
#
# Build the native macOS desktop bundle and release zip.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "============================================================"
echo "  Building Fraggler Diagnostics for macOS"
echo "============================================================"
echo "  Project root: $PROJECT_ROOT"
echo ""

cd "$PROJECT_ROOT"

# Activate venv
if [ -d "fraggler-mac310-venv" ]; then
    source fraggler-mac310-venv/bin/activate
else
    echo "ERROR: fraggler-mac310-venv not found."
    exit 1
fi

python -m pip install -r requirements.txt -r packaging/build-requirements.txt

rm -rf build dist

echo ""
echo "Running unified desktop build..."
echo ""

python build_qt.py

echo ""
echo "============================================================"
echo "  ✅ Build complete!"
echo "  App bundle : dist/Fraggler.app"
echo "  Release zip: dist/releases/Fraggler_macOS.zip"
echo ""
echo "  App size   : $(du -sh dist/Fraggler.app | cut -f1)"
echo "============================================================"
