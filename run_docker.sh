#!/bin/bash

# --- Color Definitions ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}============================================================${NC}"
echo -e "${BLUE}   Fraggler Diagnostics — Docker Runtime (GUI)${NC}"
echo -e "${BLUE}============================================================${NC}"

# 1. Build the Docker Image
echo -e "${GREEN}[1/2] Building Docker image...${NC}"
docker build -f packaging/Dockerfile.runtime -t fraggler-runtime .

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Docker build failed.${NC}"
    exit 1
fi

# 2. Allow Docker to connect to X11 (local)
echo -e "${GREEN}[2/2] Granting X11 permissions and launching app...${NC}"
xhost +local:docker > /dev/null

# 3. Run the container
# -v /tmp/.X11-unix:/tmp/.X11-unix: Map the X11 socket
# -e DISPLAY: Pass the display variable
# --device /dev/dri: Allow GPU acceleration (optional but helpful)
docker run --rm -it \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --name fraggler-app \
    fraggler-runtime

# 4. Revoke permissions on exit
xhost -local:docker > /dev/null
echo -e "${BLUE}App closed.${NC}"
