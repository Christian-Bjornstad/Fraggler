#!/bin/bash
# Download all Python dependencies as wheels for offline Linux installation
# Run this on a machine with internet

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/linux_offline_deps"
RUNTIME_REQ="$PROJECT_ROOT/requirements.txt"
BUILD_REQ="$SCRIPT_DIR/build-requirements.txt"
PYTHON_VERSION="3.10"

echo "=========================================="
echo "  Downloading wheels for offline Linux"
echo "=========================================="

mkdir -p "$OUTPUT_DIR"
cd "$PROJECT_ROOT"

echo "Downloading wheels to $OUTPUT_DIR..."
python3 -m pip download \
  -r "$RUNTIME_REQ" \
  -r "$BUILD_REQ" \
  -d "$OUTPUT_DIR" \
  --platform manylinux_2_17_x86_64 \
  --platform linux_x86_64 \
  --only-binary=:all: \
  --python-version "$PYTHON_VERSION"

echo ""
echo "Download complete!"
echo ""
echo "To install on offline Linux machine:"
echo "  1. Copy 'linux_offline_deps/' folder to Linux machine"
echo "  2. Create/activate a Python 3.10 virtual environment"
echo "  3. Run: pip install --no-index --find-links=linux_offline_deps/ -r requirements.txt"
echo ""
echo "To build on Linux:"
echo "  1. Install system dependencies (see linux_system_deps.sh)"
echo "  2. Copy source code to Linux"
echo "  3. Run: python build_qt.py"
