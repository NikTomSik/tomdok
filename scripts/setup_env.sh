#!/usr/bin/env bash
# Install/verify WRN-repurposing environment. Idempotent.
#   bash setup_env.sh            # install missing
#   bash setup_env.sh --check    # check only, no installs
set -uo pipefail
MODE="${1:-install}"
ENV=screen; PY=3.11
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniconda3}"

dl(){ local url=$1 out=$2
  if [ -x /usr/bin/curl ]; then /usr/bin/curl -fsSL "$url" -o "$out" && return 0; fi
  curl -fsSL "$url" -o "$out" && return 0
  python3 -c "import urllib.request,sys; urllib.request.urlretrieve(sys.argv[1],sys.argv[2])" "$url" "$out"; }

# 1) Miniconda
if [ "$MODE" != "--check" ] && ! command -v conda >/dev/null 2>&1 && [ ! -x "$CONDA_ROOT/bin/conda" ]; then
  echo "[setup] Installing Miniconda..."; dl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh /tmp/m.sh
  bash /tmp/m.sh -b -p "$CONDA_ROOT"; fi
[ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ] && source "$CONDA_ROOT/etc/profile.d/conda.sh"
command -v conda >/dev/null 2>&1 || { echo "conda unavailable"; exit 1; }

# 2) environment
if [ "$MODE" != "--check" ] && ! conda env list | awk '{print $1}' | grep -qx "$ENV"; then
  echo "[setup] Creating $ENV..."; conda create -y -n "$ENV" python="$PY"; fi
source "$CONDA_ROOT/etc/profile.d/conda.sh"; conda activate "$ENV"

have_py(){ python -c "import $1" 2>/dev/null; }
have_bin(){ command -v "$1" >/dev/null 2>&1; }
try_conda(){ conda install -y -c conda-forge -c bioconda "$@" >/dev/null 2>&1; }
try_pip(){ python -m pip install -q "$@" >/dev/null 2>&1; }

install_py(){
  if have_py "$1"; then echo "OK  py  $1"; return; fi
  if [ "$MODE" = "--check" ]; then echo "MISS py $1 -> conda:$2 / pip:$3"; return; fi
  try_conda "$2" && { have_py "$1" && { echo "OK  py  $1 (conda:$2)"; return; }; }
  try_pip  "$3" && { have_py "$1" && { echo "OK  py  $1 (pip:$3)";  return; }; }
  echo "MISS py $1 -> conda:$2 / pip:$3"; }
install_bin(){
  if have_bin "$1"; then echo "OK  bin $1"; return; fi
  if [ "$MODE" = "--check" ]; then echo "MISS bin $1 -> conda:$2"; return; fi
  try_conda "$2" && { have_bin "$1" && { echo "OK  bin $1 (conda:$2)"; return; }; }
  [ -n "${3:-}" ] && { try_pip "$3" && { have_bin "$1" && { echo "OK  bin $1 (pip:$3)"; return; }; }; }
  echo "MISS bin $1 -> conda:$2"; }

echo "[setup] Mode: $MODE"
echo "[setup] Python dependencies..."
install_py rdkit rdkit rdkit
install_py numpy numpy numpy
install_py pdbfixer pdbfixer pdbfixer
install_py openmm openmm openmm
install_py meeko meeko meeko
install_py dimorphite_dl dimorphite-dl dimorphite-dl
install_py matplotlib matplotlib matplotlib
install_py PIL pillow Pillow
install_py plip.structure.preparation plip plip
echo "[setup] Binaries..."
install_bin obabel openbabel
install_bin gnina gnina
install_bin vina autodock-vina

# LeDock
if have_bin ledock && have_bin lepro; then echo "OK  bin ledock/lepro";
elif [ "$MODE" != "--check" ]; then
  echo "[setup] Attempting LeDock install..."
  dl "http://www.lephar.com/software/ledock_linux_x86_64.tar.gz" /tmp/ledock.tgz 2>/dev/null \
    && { cd /tmp && tar -xzf ledock.tgz 2>/dev/null \
         && install -m755 /tmp/ledock "$CONDA_ROOT/envs/$ENV/bin/ledock" 2>/dev/null || true \
         && install -m755 /tmp/lepro  "$CONDA_ROOT/envs/$ENV/bin/lepro" 2>/dev/null || true; }
  have_bin ledock && have_bin lepro && echo "OK  bin ledock/lepro" || echo "MISS bin ledock/lepro -> manually from http://www.lephar.com/download.htm"
else echo "MISS bin ledock/lepro (check)"; fi

echo; echo "=== SUMMARY ($MODE) ==="
echo "[setup] Done. Activate: conda activate $ENV"
