#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
HOST_UID="$(id -u)"
HOST_GID="$(id -g)"

echo "============================================================"
echo "  Building Fraggler Diagnostics for Linux Offline Bundle"
echo "============================================================"

# Ensure Docker binaries are in PATH (common for Docker Desktop on Mac)
export PATH="/Applications/Docker.app/Contents/Resources/bin:/usr/local/bin:$PATH"

if ! command -v docker &> /dev/null; then
    echo "ERROR: docker command not found in PATH."
    exit 1
fi

cd "$PROJECT_ROOT"

echo "Building Docker image..."
docker build -f packaging/Dockerfile.linux -t fraggler-linux-build .

echo "Running build..."
mkdir -p "$PROJECT_ROOT/dist"
docker run --rm \
    -v "$PROJECT_ROOT/dist:/mnt" \
    alpine sh -lc "rm -rf /mnt/* /mnt/.[!.]* /mnt/..?* 2>/dev/null || true"
docker run --rm \
    -v "$PROJECT_ROOT/dist:/app/dist" \
    fraggler-linux-build \
    python3 build_qt.py
docker run --rm \
    -v "$PROJECT_ROOT/dist:/mnt" \
    alpine sh -lc "chown -R ${HOST_UID}:${HOST_GID} /mnt"

echo "Done!"
echo "Portable folder: dist/Fraggler_Linux"
echo "Offline zip    : dist/releases/Fraggler_Linux_offline.zip"
