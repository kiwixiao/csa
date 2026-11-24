#!/bin/bash
# Smoke test of every tool the chain uses; prints versions and OK/FAIL per item. Exit 1 if anything fails.
#   bash cfd_postproc/check_env.sh        (after: source cfd_postproc/activate.sh, or on the original machine as is)
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(dirname "$HERE"); FAIL=0
chk() { local name=$1; shift; local out; if out=$("$@" 2>&1); then echo "OK    $name  $(echo "$out" | tail -n 1)"; else echo "FAIL  $name"; echo "$out" | tail -n 3 | sed 's/^/      /'; FAIL=1; fi; }
chk "pvbatch (ParaView)" bash -c "'$HERE/pvb.sh' -c \"import paraview, numpy, scipy; from paraview.simple import OpenDataFile; print('paraview', paraview.__version__, 'numpy', numpy.__version__)\" 2>/dev/null | grep paraview"
chk "python3 (plots)" python3 -c "import sys, numpy, scipy, matplotlib, pandas, PIL; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'matplotlib', matplotlib.__version__, 'pandas', pandas.__version__)"
chk "ffmpeg" bash -c "ffmpeg -version | head -1 && ffmpeg -hide_banner -encoders 2>/dev/null | grep -q libx264"
chk "node + docx" bash -c "cd '$HERE' && node -e \"require('docx'); console.log('node', process.version, 'docx', JSON.parse(require('fs').readFileSync('node_modules/docx/package.json')).version)\""
chk "mirtk (env mirtk)" bash -c "conda run -n mirtk mirtk help >/dev/null && echo mirtk ok"
chk "csa-bifurcation (env mirtk)" bash -c "conda run -n mirtk csa-bifurcation --help >/dev/null && conda run -n mirtk pip show csa-slicer | grep -i location"
chk "vmtk (env vmtk)" bash -c "conda run -n vmtk python -c \"from vmtk import vmtkscripts; import vtk; print('vmtk ok, vtk', vtk.VTK_VERSION)\""
[ $FAIL = 0 ] && echo "all OK" || echo "SOME CHECKS FAILED"
exit $FAIL
