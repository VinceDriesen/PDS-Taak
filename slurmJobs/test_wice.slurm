#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --account=lp_h_pds_iiw
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=00:20:00
#SBATCH --output="%x_%j.out"
#SBATCH --error="%x_%j.err"

source /data/leuven/303/vsc30380/slurmhooks

module load CMake
module load CUDA

REPO_URL="https://github.com/VinceDriesen/PDS-Taak.git"
REPO_ROOT="my_gpu_project"    # The folder git creates
SUB_DIR="ParticleSim"         # The folder containing CMakeLists.txt
BUILD_DIR="build"             # The build directory
EXE_PATH="bin/particlesim"    # The binary location (relative to build dir)

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

# 3. Enter subdirectory and Build
if [ -d "$SUB_DIR" ]; then
    cd "$SUB_DIR"
else
    echo "Error: Directory $SUB_DIR not found!"
    exit 1
fi

echo "Configuring CMake..."
cmake -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release .

echo "Compiling..."
cmake --build "$BUILD_DIR"

echo "Starting simulation..."
./"$EXE_PATH" test nocpu noio