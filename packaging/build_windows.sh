#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "============================================================"
echo "  Building Fraggler Diagnostics for Windows (.exe via Docker)"
echo "============================================================"

# Ensure Docker binaries are in PATH (common for Docker Desktop on Mac)
export PATH="/Applications/Docker.app/Contents/Resources/bin:/usr/local/bin:$PATH"

if ! command -v docker &> /dev/null; then
    echo "ERROR: docker command not found in PATH."
    exit 1
fi

cd "$PROJECT_ROOT"

echo "Building Docker image..."
docker build -f packaging/Dockerfile.windows -t fraggler-windows-build .

echo "Running build..."
# We map the dist/ folder to export the compiled executable
docker run --rm -v "$PROJECT_ROOT/dist:/app/dist" fraggler-windows-build

echo "Done! Windows executable is available in dist/"
