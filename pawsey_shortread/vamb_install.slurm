#!/bin/bash
#SBATCH --job-name=vamb_setup
#SBATCH --partition=gpu-dev
#SBATCH --gres=gpu:1
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# submit this with sbatch --account=${PAWSEY_PROJECT}-gpu vamb_install.slurm

###### NOTES ######
#
# This script creates a virtual environment for vamb on a GPU node at Pawsey.
# There are two (or 3?) things that you need to get working together.
#
# vamb needs a specific version of pytorch, and so you need to make
# sure that the rocm build has that. 
#
# You can check the pytorch versions at https://download.pytorch.org/whl/rocm6.2.4 
# but I found it a bit of trial and error to get the right combination.
#
# First, load the right module.
# Second, install the right pytorch version.
# Third, we install vamb, but with no dependencies, otherwise it will remove our rocm pytorch
# and install the NVIDIA one.
# Fourth, we install the dependencies for vamb, except for torch.
# Finally, we use pip check to list any missing dependencies.
#
# When vamb gets updated, you need to look for the torch version that include +rocm
# at the end. For example, were using torch-2.6.0+rocm6.2.4-cp312 because we have
# torch 2.6.0, rocm 6.2.4 python 3.12.
#
#



set -euo pipefail

# ---- 1) Get the right module here!
module load rocm/6.2.4

# ---- 2) Create and load the virtual environment 
VENV_DIR="/scratch/pawsey1018/edwa0468/software/pip/vamb"
if [[ -e "$VENV_DIR" ]]; then 
	rm -rf "$VENV_DIR"
fi

mkdir -p "$VENV_DIR"
python -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# Upgrade packaging tooling
pip install --upgrade pip setuptools wheel

# ---- 3) Install pytorch and then vamb
# The pytorch module already provides torch/ROCm; avoid pip installing another torch build.
# This requires pip > 23.2
pip --version
## pip install vamb --exclude torch
pip install --index-url https://download.pytorch.org/whl/rocm6.2.4 "torch == 2.6.0" torchvision torchaudio
pip install --no-deps  vamb "vambcore == 0.1.2" "numpy == 1.26.4" "pycoverm == 0.6.2" "networkx == 3.4.2" "scikit-learn == 1.6.1" "dadaptation == 3.2" "loguru == 0.7.3" "pyhmmer == 0.10.15" "pyrodigal == 3.6.3" joblib scipy threadpoolctl archspec psutil~=5.8

echo "Running pip check"
pip check

echo "Sane?"
# ---- Sanity checks ----
python - <<'PY'
import sys
print("Python:", sys.version)
try:
    import torch
    print("PyTorch:", torch.__version__)
    print("CUDA available? (ROCm builds return True here):", torch.cuda.is_available())
    if torch.cuda.is_available():
        print("Device 0:", torch.cuda.get_device_name(0))
except Exception as e:
    print("Torch check failed:", e)

try:
    import vamb
    print("VAMB:", getattr(vamb, "__version__", "unknown"))
except Exception as e:
    print("VAMB import failed:", e)
PY

echo "Done. Virtual env at: $VENV_DIR"

