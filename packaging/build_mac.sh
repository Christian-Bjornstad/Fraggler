#!/bin/bash
# ============================================================
# Build Fraggler Diagnostics — Mac Executable
# ============================================================
# Prerequisites:
#   - macOS with the fraggler-mac310-venv virtualenv
#   - Run from the OUS/ project root
#
# Usage:
#   ./packaging/build_mac.sh
#
# Output:
#   packaging/dist/fraggler-diagnostics/
# ============================================================

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

# Install PyInstaller if not present (prefer 6.x)
pip install pyinstaller 2>/dev/null || true

# Clean previous builds
rm -rf packaging/build packaging/dist

# Build using the dis.py monkey-patch wrapper
# (Workaround for Python 3.10.0 dis module bug)
echo ""
echo "Running PyInstaller..."
echo ""

python -c "
import dis
_orig = dis._get_const_info
def _patched(ci, cl):
    try:
        return _orig(ci, cl)
    except IndexError:
        return ci, repr(ci)
dis._get_const_info = _patched

import PyInstaller.__main__
PyInstaller.__main__.run([
    'packaging/fraggler_diagnostics.spec',
    '--distpath', 'packaging/dist',
    '--workpath', 'packaging/build',
    '--clean',
    '--noconfirm',
])
"

echo ""
echo "============================================================"
echo "  ✅ Build complete!"
echo "  Executable: packaging/dist/fraggler-diagnostics/"
echo ""
echo "  Size: $(du -sh packaging/dist/fraggler-diagnostics/ | cut -f1)"
echo ""
echo "  To run:"
echo "    ./packaging/dist/fraggler-diagnostics/fraggler-diagnostics"
echo "============================================================"
