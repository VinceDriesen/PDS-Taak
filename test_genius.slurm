#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --account=lp_h_pds_iiw
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:20:00
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

# 1. Load Environment
source /data/leuven/303/vsc30380/slurmhooks

# CUDA 11.7 works well for both P100 (Pascal) and V100 (Volta)
module load CMake
module load CUDA

# Check versions
nvcc --version
cmake --version

# --- CONFIGURATION ---
REPO_URL="https://github.com/VinceDriesen/PDS-Taak.git"
REPO_ROOT="my_gpu_project"    # The folder git creates
SUB_DIR="ParticleSim"         # The folder containing CMakeLists.txt
BUILD_DIR="build"             # The build directory
EXE_PATH="bin/particlesim"    # The binary location (relative to build dir)
# ---------------------

# 2. Clone or Update
if [ -d "$REPO_ROOT" ]; then
    echo "Updating repo..."
    cd "$REPO_ROOT"
    git pull
else
    echo "Cloning repo..."
    git clone "$REPO_URL" "$REPO_ROOT"
    cd "$REPO_ROOT"
fi

# 3. Build
# 3. Enter subdirectory and Build
if [ -d "$SUB_DIR" ]; then
    cd "$SUB_DIR"
else
    echo "Error: Directory $SUB_DIR not found!"
    exit 1
fi

echo "Configuring CMake..."
# Generate build files
cmake -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release .

echo "Compiling..."
# Build using 18 CPUs for speed
cmake --build "$BUILD_DIR" --parallel "$SLURM_CPUS_PER_TASK"

# 4. Run
echo "Starting simulation..."
# The executable is inside build/bin/particlesim based on your description
./"$SUB_DIR"/"$EXE_PATH" test nocpi noio