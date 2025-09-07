#!/bin/bash
#SBATCH --job-name=vamb_setup
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# submit this with sbatch --account=${PAWSEY_PROJECT}-gpu vamb_install.slurm

set -euo pipefail

# ---- 1) Load GPU/ROCm module ----
# module load rocm/6.4.1 <- note that this doesn't have the correct version of pytorch required by vamb
module load rocm/6.2.4

# (Optional but often helpful)
echo "=== Module list ==="
module list || true

# ---- 2) Create virtual environment ----
# NOTE: We need pip >23.2. At the time of writing, the default python module is 3.11.4, which has pip 23.0.1.
# so we make the venv separately.
VENV_DIR="/scratch/pawsey1018/edwa0468/software/pip/vamb"
if [[ -e "$VENV_DIR" ]]; then 
	rm -rf "$VENV_DIR"
fi

mkdir -p "$VENV_DIR"
python -m venv "$VENV_DIR"


# Activate venv
source "$VENV_DIR/bin/activate"

# Upgrade packaging tooling
pip install --upgrade pip setuptools wheel

# ---- 3) Install vamb (without dragging a conflicting torch) ----
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

