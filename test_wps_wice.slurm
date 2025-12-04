#!/bin/bash
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

APP = ./"$EXE_PATH" noio gpu

# ---------------------------------------------------------
# 2. MPS Configuratie & Opstarten
# ---------------------------------------------------------

# Maak unieke tijdelijke directories aan voor deze job (voorkomt conflicten)
export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-pipe-$SLURM_JOB_ID
export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-mps-log-$SLURM_JOB_ID
mkdir -p $CUDA_MPS_PIPE_DIRECTORY
mkdir -p $CUDA_MPS_LOG_DIRECTORY

echo "Starten van NVIDIA MPS Daemon..."
# Start de MPS control daemon
nvidia-cuda-mps-control -d

# Wacht even om zeker te zijn dat de daemon draait
sleep 2

# Controleer of MPS actief is
if pgrep "nvidia-cuda-mps" > /dev/null
then
    echo "MPS Daemon is succesvol gestart."
else
    echo "FOUT: MPS Daemon kon niet starten."
    exit 1
fi

# ---------------------------------------------------------
# 3. De Loop: 20, 40, 60, 80, 100
# ---------------------------------------------------------

for PCT in 20 40 60 80 100
do
    echo "=========================================="
    echo "RUN MET GPU LIMIET: $PCT%"
    echo "=========================================="

    # De cruciale regel uit jouw strategie:
    export CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=$PCT

    # Voer de applicatie uit (gebruik 'time' om prestaties te meten)
    # Zorg dat je applicatie uitvoer geeft die je later kunt analyseren
    /usr/bin/time -f "Real time: %E sec" $APP

    # Korte pauze tussen runs (optioneel, voor stabiliteit)
    sleep 2
done

# ---------------------------------------------------------
# 4. Opruimen (Teardown)
# ---------------------------------------------------------

echo "Stoppen van MPS Daemon..."
echo quit | nvidia-cuda-mps-control

# Verwijder tijdelijke mappen
rm -rf $CUDA_MPS_PIPE_DIRECTORY
rm -rf $CUDA_MPS_LOG_DIRECTORY

echo "Klaar."