#!/bin/bash
# One-command setup of everything the CSA + CFD post-processing chain needs (Linux x86_64, conda >= 23 on PATH).
# Idempotent: anything already present is left as it is.
#
#   bash cfd_postproc/setup.sh [--paraview-prefix DIR] [--mirtk-install DIR]
#
#   1. conda env 'mirtk'   MIRTK + mirtk-pipeline (mirtk-install/install.sh), then csa-slicer from THIS repo (editable)
#   2. conda env 'vmtk'    vmtk 1.5.0 (env/vmtk_env.yml); csa-bifurcation calls it by this name
#   3. conda env 'cfdpost' python 3.11 + numpy/scipy/matplotlib/pandas/pillow + ffmpeg + node (env/cfdpost_env.yml)
#   4. ParaView 6.0.1 (binary download into --paraview-prefix, default ~/opt) -> pvbatch for pvb.sh
#   5. npm install (docx) in cfd_postproc
# Afterwards:  source cfd_postproc/activate.sh   then   bash cfd_postproc/check_env.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(dirname "$HERE")
PV_PREFIX=$HOME/opt; MIRTK_INSTALL=""
while [ $# -gt 0 ]; do
  case $1 in
    --paraview-prefix) PV_PREFIX=$2; shift 2 ;;
    --mirtk-install) MIRTK_INSTALL=$2; shift 2 ;;
    *) echo "unknown option $1"; exit 1 ;;
  esac
done
command -v conda >/dev/null || { echo "conda not found: install Miniconda first"; exit 1; }
has_env() { conda env list | awk '{print $1}' | grep -qx "$1"; }
say() { echo; echo "=== $* ==="; }

say "1. conda env mirtk (MIRTK + registration pipeline + CSA slicer)"
if has_env mirtk; then
  echo "env mirtk exists - left as is"
else
  if [ -z "$MIRTK_INSTALL" ]; then
    MIRTK_INSTALL=$(dirname "$REPO")/mirtk-install
    [ -d "$MIRTK_INSTALL" ] || git clone https://github.com/kiwixiao/mirtk-install.git "$MIRTK_INSTALL"
  fi
  bash "$MIRTK_INSTALL/install.sh"
fi
LOC=$(conda run -n mirtk pip show csa-slicer 2>/dev/null | awk '/^Editable project location/ {print $4}')
if [ "$LOC" != "$REPO/python_slicer" ]; then
  echo "csa-slicer -> editable install from $REPO/python_slicer (was: ${LOC:-not editable / GitHub copy})"
  conda run -n mirtk pip install -e "$REPO/python_slicer"
else
  echo "csa-slicer already runs from $LOC"
fi

say "2. conda env vmtk"
has_env vmtk && echo "env vmtk exists - left as is" || conda env create -f "$HERE/env/vmtk_env.yml"

say "3. conda env cfdpost (plots / tables / videos)"
has_env cfdpost && echo "env cfdpost exists - left as is" || conda env create -f "$HERE/env/cfdpost_env.yml"

say "4. ParaView 6.0.1"
PV=ParaView-6.0.1-MPI-Linux-Python3.12-x86_64
if [ -x /opt/$PV/bin/pvbatch ]; then
  PVBATCH=/opt/$PV/bin/pvbatch
elif [ -x "$PV_PREFIX/$PV/bin/pvbatch" ]; then
  PVBATCH=$PV_PREFIX/$PV/bin/pvbatch
else
  mkdir -p "$PV_PREFIX"
  curl -L -o "$PV_PREFIX/$PV.tar.gz" "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v6.0&type=binary&os=Linux&downloadFile=$PV.tar.gz"
  tar -xzf "$PV_PREFIX/$PV.tar.gz" -C "$PV_PREFIX" && rm "$PV_PREFIX/$PV.tar.gz"
  PVBATCH=$PV_PREFIX/$PV/bin/pvbatch
fi
echo "pvbatch: $PVBATCH"

say "5. node package (docx)"
( cd "$HERE" && conda run --no-capture-output -n cfdpost npm install --no-fund --no-audit )

cat > "$HERE/activate.sh" <<EOT
# source this before running the pipeline (written by setup.sh)
eval "\$(conda shell.bash hook)"
conda activate cfdpost
export PVBATCH=$PVBATCH
EOT
echo; echo "Done. Next:  source $HERE/activate.sh && bash $HERE/check_env.sh"
