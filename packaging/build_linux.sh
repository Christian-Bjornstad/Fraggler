#!/bin/bash
# ============================================================
# Build Fraggler Diagnostics — Linux Executable (via Docker)
# ============================================================
# Prerequisites:
#   - Docker installed and running
#   - Run from the OUS/ project root
#
# Usage:
#   ./packaging/build_linux.sh
#
# Output:
#   packaging/dist/fraggler-diagnostics-linux/
# ============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "============================================================"
echo "  Building Fraggler Diagnostics for Linux (via Docker)"
echo "============================================================"
echo "  Project root: $PROJECT_ROOT"
echo ""

cd "$PROJECT_ROOT"

# Ensure Docker binaries are in PATH (common for Docker Desktop on Mac)
export PATH="/Applications/Docker.app/Contents/Resources/bin:/usr/local/bin:$PATH"

if ! command -v docker &> /dev/null; then
    echo "ERROR: docker command not found in PATH."
    exit 1
fi

if ! docker info &> /dev/null; then
    echo "ERROR: Docker daemon is not running."
    echo "Please ensure Docker Desktop is open and the whale icon is visible."
    exit 1
fi
echo "Using Docker at: $(command -v docker)"

# Create output dir
mkdir -p packaging/dist

# Build Docker image
echo "Building Docker image..."
docker build \
    -f packaging/Dockerfile.linux \
    -t fraggler-diagnostics-build \
    .

# Run the build inside Docker and copy output
echo ""
echo "Running PyInstaller inside Docker..."
docker run --rm \
    -v "$PROJECT_ROOT/packaging/dist:/output" \
    fraggler-diagnostics-build

echo ""
echo "============================================================"
echo "  ✅ Build complete!"
echo "  Executable: packaging/dist/fraggler-diagnostics-linux/"
echo ""
echo "  To test with Docker:"
echo "    docker run --rm -it -p 5078:5078 \\"
echo "      -v packaging/dist/fraggler-diagnostics-linux:/app \\"
echo "      python:3.10-slim /app/fraggler-diagnostics"
echo "============================================================"
