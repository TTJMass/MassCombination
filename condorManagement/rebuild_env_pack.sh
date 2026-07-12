#!/usr/bin/env bash
# Rebuilds envCache/masscomb_packed.tar.gz, the conda-pack tarball shipped to
# every Condor job (see createCondorJobs.py --conda-pack-tarball) so jobs
# don't have to activate conda directly from EOS.
#
# Run this whenever pyconvino source changes, or whenever a package in the
# masscomb conda env is added/updated/removed. The packed tarball is a point
# -in-time snapshot; createCondorJobs.py does not rebuild it automatically.
#
# conda-pack refuses to pack an env with editable-installed packages, and
# convino_jax (pyconvino) is normally installed editable for interactive dev.
# This script temporarily reinstalls it non-editable, packs, then restores
# the editable install (via a trap, so it runs even if packing fails).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYCONVINO_DIR="$REPO_ROOT/pyconvino"
ENV_CACHE_DIR="$REPO_ROOT/envCache"
PACKED_TARBALL="$ENV_CACHE_DIR/masscomb_packed.tar.gz"
CONDA_SH="/eos/home-s/sewuchte/MyConda/etc/profile.d/conda.sh"

if [ -f "$CONDA_SH" ]; then
    source "$CONDA_SH"
else
    export PATH="/eos/home-s/sewuchte/MyConda/bin:$PATH"
fi
conda activate masscomb
export PYTHONNOUSERSITE=1

mkdir -p "$ENV_CACHE_DIR"

restore_editable() {
    echo "Restoring editable convino_jax install..."
    pip uninstall -y convino_jax >/dev/null 2>&1 || true
    pip install -e "$PYCONVINO_DIR" --no-deps --no-user
}
trap restore_editable EXIT

echo "Ensuring conda-pack is installed..."
python3 -c "import conda_pack" 2>/dev/null || pip install --no-user conda-pack

echo "Temporarily installing convino_jax non-editable (required for conda-pack)..."
pip install "$PYCONVINO_DIR" --no-deps --no-user --force-reinstall

echo "Packing masscomb env to $PACKED_TARBALL ..."
conda-pack -n masscomb -o "$PACKED_TARBALL" --force

echo "Done."
ls -lh "$PACKED_TARBALL"
