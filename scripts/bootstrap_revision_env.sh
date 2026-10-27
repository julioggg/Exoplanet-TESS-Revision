#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${1:-$HOME/exoplanet_revision/repo}"
VENV_DIR="${EXOPLANET_REVISION_VENV:-$HOME/venvs/exoplanet-revision}"

cd "$REPO_DIR"

echo "============================================================"
echo "EXOPLANET REVISION ENVIRONMENT BOOTSTRAP"
echo "============================================================"
echo "Repository : $REPO_DIR"
echo "Virtualenv : $VENV_DIR"
echo "Python     : $(command -v python)"
python --version
echo

mkdir -p "$(dirname "$VENV_DIR")"

if [ ! -x "$VENV_DIR/bin/python" ]; then
    echo "[1/4] Creating virtual environment..."
    python -m venv "$VENV_DIR"
else
    echo "[1/4] Virtual environment already exists."
fi

echo "[2/4] Upgrading packaging tools..."
"$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel

echo "[3/4] Installing revision requirements..."
"$VENV_DIR/bin/python" -m pip install -r requirements_revision.txt

echo "[4/4] Running dependency check..."
"$VENV_DIR/bin/python" revision/99_check_environment.py

"$VENV_DIR/bin/python" -m pip freeze > environment_revision_freeze.txt

echo
echo "Exact package versions written to:"
echo "  $REPO_DIR/environment_revision_freeze.txt"
echo
echo "Activate with:"
echo "  source $VENV_DIR/bin/activate"
echo "============================================================"
