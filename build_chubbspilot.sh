#!/bin/bash

set -x

# Function to determine the local network range
get_network_range() {
    local ip
    ip=$(ipconfig getifaddr en0 2>/dev/null || ip addr show | grep 'inet ' | grep -v '127.0.0.1' | awk '{print $2}' | cut -d/ -f1 | head -n1)
    if [[ -z "$ip" ]]; then
        echo "ERROR: Unable to determine local IP address."
        exit 1
    fi
    # Convert IP to subnet (e.g., 192.168.1.0/24)
    local subnet
    subnet=$(echo "$ip" | awk -F. '{print $1 "." $2 "." $3 ".0/24"}')
    echo "$subnet"
}

# Function to detect device IP using nmap
get_device_ip_nmap() {
    local network_range
    network_range=$(get_network_range)
    echo "Scanning network range: $network_range"
    # Scan for devices with SSH ports 22 or 8022 open
    nmap -p 22,8022 --open -oG - "$network_range" | awk '/22\/open|8022\/open/{print $2}'
}

# Automatically detect DEVICE_IP
detect_device_ip() {
    local possible_ips
    possible_ips=$(get_device_ip_nmap)
    local device_ip=""
    for ip in $possible_ips; do
        if ssh -o ConnectTimeout=2 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "$DEVICE_USER@$ip" 'exit' 2>/dev/null; then
            device_ip="$ip"
            break
        fi
    done

    if [[ -z "$device_ip" ]]; then
        echo "ERROR: Failed to detect device IP via nmap."
        exit 1
    fi

    echo "$device_ip"
}

# Configuration
DEVICE_USER="comma"
DEVICE_IP=$(detect_device_ip)

# Remove existing SSH keys for DEVICE_IP
ssh-keygen -R "$DEVICE_IP" > /dev/null 2>&1

echo "Detected Device IP: $DEVICE_IP"

# Function to get the current git branch on the SSH device
get_device_branch() {
    ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "$DEVICE_USER@$DEVICE_IP" 'cd /data/openpilot && git rev-parse --abbrev-ref HEAD'
}

# Function to check if prebuilt exists
check_prebuilt() {
    # Check for the 'prebuilt' file
    PREBUILT_FILE=$(ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "$DEVICE_USER@$DEVICE_IP" 'test -f /data/openpilot/prebuilt && echo "yes" || echo "no"')

    # Additionally, check if any .so files exist as an indicator
    ESSENTIAL_SO=$(ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "$DEVICE_USER@$DEVICE_IP" 'find /data/openpilot -type f -name "*.so" | grep -q . && echo "yes" || echo "no"')

    if [[ "$PREBUILT_FILE" == "yes" || "$ESSENTIAL_SO" == "yes" ]]; then
        echo "yes"
    else
        echo "no"
    fi
}

# Function to skip compilation and clean
skip_compile_device() {
    echo "Proceeding without recompilation..." >&2  # Redirect to stderr
    OUTPUT=$(ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null "$DEVICE_USER@$DEVICE_IP" << 'ENDSSH'
        set -e
        set -x
        cd /data/openpilot

        # Add compiled files
        git add -f .
        git commit -m "compile" --allow-empty

        # Delete all .cc files except specified ones
        find . -type f -name '*.cc' ! \( -name 'main.cc' -o -name 'map_renderer.cc' -o -name 'transform.cc' \) -delete

        # Remove other unwanted file types
        find . -name '*.o' -delete
        find . -name '*.os' -delete
        find . -name '*.pyc' -delete
        find . -name 'moc_*' -delete
        find . -name '__pycache__' -delete
        find . -name '*.a' -delete

        # Remove build system artifacts
        rm -rf .sconsign.dblite Jenkinsfile release/

        # Remove generated directories
        rm -rf selfdrive/ubloxd/generated cereal/gen/cpp

        # Remove remaining MPC-generated files
        rm -rf selfdrive/controls/lib/longitudinal_mpc_lib/c_generated_code/main_*
        rm -rf selfdrive/controls/lib/lateral_mpc_lib/c_generated_code/main_*

        # Additional Deletions
        # Remove specific directories except .github/workflows and delete body/
        find .github -mindepth 1 ! -path './.github/workflows*' -delete
        rm -rf .vscode .devcontainer debug \
               tools/cabana tools/camerastream tools/car_porting tools/joystick \
               tools/latencylogger tools/lib/tests tools/plotjuggler tools/profiling \
               tools/replay tools/scripts tools/serial tools/sim tools/ssh \
               tools/ubuntu_setup tools/webcam body/

        # Remove Docker-related files
        rm -f Dockerfile* .dockerignore

        # Remove test files and directories
        find . -type f -iname '*test*' -delete
        find . -type d -iname '*test*' -exec rm -rf {} +

        # Remove all files named LICENSE (case-insensitive)
        find . -type f -iname 'license*' -delete

        # Remove teleoprtc files
        find . -type f -iname '*teleoprtc*' -delete

        # Remove configuration files
        rm -f .clang-tidy .editorconfig

        # Remove debug-related files and directories
        find . -type f -iname '*debug*' -delete
        find . -type d -iname '*debug*' -exec rm -rf {} +

        touch prebuilt

        # Amend with cleaned state
        git add -f .
        git commit --amend -m "Prebuilt" --allow-empty

        # Output hash
        echo "HASH=$(git rev-parse HEAD)"
ENDSSH
    )
    set -x  # Re-enable debug mode

    SSH_STATUS=$?
    if [ $SSH_STATUS -ne 0 ]; then
        echo "ERROR: Commit and clean SSH session failed with status $SSH_STATUS" >&2
        echo "Output was:" >&2
        echo "$OUTPUT" >&2
        exit 1
    fi

    # Extract hash
    COMMIT_HASH=$(echo "$OUTPUT" | grep "^HASH=" | cut -d= -f2)

    if ! [[ $COMMIT_HASH =~ ^[0-9a-f]{40}$ ]]; then
        echo "ERROR: Invalid commit hash: $COMMIT_HASH" >&2
        echo "Full output was:" >&2
        echo "$OUTPUT" >&2
        exit 1
    fi

    echo "$COMMIT_HASH"
}

# Function to handle workspace operations
workspace_operations() {
    local device_commit=$1
    echo "=== Starting workspace operations ==="
    export GIT_SSH_COMMAND="ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"

    # Reset workspace completely
    git reset --hard HEAD
    git clean -fd

    # Setup remote if needed
    if (! git remote | grep -q "device"); then
        git remote add device ssh://"${DEVICE_USER}"@"${DEVICE_IP}":/data/openpilot
    fi

    # Fetch and checkout
    git fetch device "$BRANCH_NAME" --force
    git checkout "$BRANCH_NAME" || { echo "ERROR: Failed to checkout branch '$BRANCH_NAME'"; exit 1; }

    # Cherry pick with theirs strategy
    git cherry-pick --no-commit --strategy=recursive -X theirs "$device_commit" || { echo "ERROR: Cherry-pick failed"; exit 1; }

    # Create initial prebuilt commit
    git add -f .
    git commit -m "Prebuilt $(date '+%Y-%m-%d')" || { echo "ERROR: Initial commit failed"; exit 1; }

    # Setup Git LFS
    git lfs uninstall
    git lfs install

    # Create .gitattributes
    cat > .gitattributes << 'EOF'
* text=auto

# Binary files
*.bin filter=lfs diff=lfs merge=lfs -text
*.o filter=lfs diff=lfs merge=lfs -text
*.so filter=lfs diff=lfs merge=lfs -text
*.so.* filter=lfs diff=lfs merge=lfs -text
*.a filter=lfs diff=lfs merge=lfs -text
*.dll filter=lfs diff=lfs merge=lfs -text
*.dylib filter=lfs diff=lfs merge=lfs -text
*.pyc filter=lfs diff=lfs merge=lfs -text

# Media files
*.onnx filter=lfs diff=lfs merge=lfs -text
*.svg filter=lfs diff=lfs merge=lfs -text
*.png filter=lfs diff=lfs merge=lfs -text
*.gif filter=lfs diff=lfs merge=lfs -text
*.ttf filter=lfs diff=lfs merge=lfs -text
*.wav filter=lfs diff=lfs merge=lfs -text

# Project specific
selfdrive/car/tests/test_models_segs.txt filter=lfs diff=lfs merge=lfs -text
system/hardware/tici/updater filter=lfs diff=lfs merge=lfs -text
selfdrive/ui/qt/spinner_larch64 filter=lfs diff=lfs merge=lfs -text
selfdrive/ui/qt/text_larch64 filter=lfs diff=lfs merge=lfs -text
selfdrive/ui/ui filter=lfs diff=lfs merge=lfs -text
EOF

    # Setup LFS tracking and migrate existing files
    git add .gitattributes
    git lfs track "."

    # Amend the commit with LFS-tracked files
    git add -f .
    git commit --amend --no-edit || { echo "ERROR: Amending commit with LFS files failed"; exit 1; }

    # Push final result
    git push --force origin "$BRANCH_NAME" || { echo "ERROR: Force push to origin failed"; exit 1; }
}

# Main execution
echo "Starting automated sync process..."

PREBUILT=$(check_prebuilt)

COMMIT_HASH=$(skip_compile_device) || exit 1

# Set BRANCH_NAME based on the device's current branch
BRANCH_NAME=$(get_device_branch)
if [[ -z "$BRANCH_NAME" ]]; then
    echo "ERROR: Failed to retrieve the device's current branch."
    exit 1
fi

workspace_operations "$COMMIT_HASH"
echo "Sync process completed."