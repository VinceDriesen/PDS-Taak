#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --account=lp_h_pds_iiw
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:20:00
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

# 1. Load Environment
source /data/leuven/303/vsc30380/slurmhooks

# Load the compilers you confirmed work
module load CMake
module load CUDA

# Load CMake (Ensure this version is available and compatible)
module load CMake

# Check versions for debugging
nvcc --version
cmake --version

# --- CONFIGURATION ---
REPO_URL="https://github.com/VinceDriesen/PDS-Taak.git"
REPO_ROOT="my_gpu_project"    # The folder git creates
SUB_DIR="ParticleSim"         # The folder containing CMakeLists.txt
BUILD_DIR="build"             # The build directory
EXE_PATH="bin/particlesim"    # The binary location (relative to build dir)
# ---------------------

# 2. Clone/Update
if [ -d "$REPO_ROOT" ]; then
    echo "Updating repository..."
    cd "$REPO_ROOT"
    git pull
else
    echo "Cloning repository..."
    git clone "$REPO_URL" "$REPO_ROOT"
    cd "$REPO_ROOT"
fi

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
./"$REPO_ROOT"/"$SUB_DIR"/"$EXE_PATH" test nocpi noio