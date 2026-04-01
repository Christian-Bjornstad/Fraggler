#!/bin/bash
# Download Linux wheels from inside a Linux Docker container.
# This is the safest option on macOS when packages like kiwisolver fail to
# resolve correctly with host-side pip download --platform flags.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/linux_offline_deps"
RUNTIME_REQ="$PROJECT_ROOT/requirements.txt"
BUILD_REQ="$SCRIPT_DIR/build-requirements.txt"
PYTHON_IMAGE="${PYTHON_IMAGE:-python:3.10-slim-bullseye}"

echo "=============================================="
echo "  Downloading Linux wheels via Docker"
echo "=============================================="
echo "Image       : $PYTHON_IMAGE"
echo "Output dir  : $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  -v "$PROJECT_ROOT:/workspace" \
  -w /workspace \
  "$PYTHON_IMAGE" \
  /bin/bash -lc "
    python -m pip install --upgrade pip &&
    mkdir -p /workspace/packaging/linux_offline_deps &&
    python -m pip download \
      -r /workspace/requirements.txt \
      -r /workspace/packaging/build-requirements.txt \
      -d /workspace/packaging/linux_offline_deps \
      --only-binary=:all:
  "

echo ""
echo "Download complete."
echo "Copy this folder to the offline Linux machine:"
echo "  packaging/linux_offline_deps/"
