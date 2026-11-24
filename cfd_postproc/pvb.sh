#!/bin/bash
# ParaView pvbatch with the anaconda/conda python env stripped (it breaks ParaView's own Python).
# Set PVBATCH to your pvbatch if ParaView 6.0.1 is installed elsewhere (see README.md, setup.sh).
unset PYTHONPATH PYTHONHOME
exec "${PVBATCH:-/opt/ParaView-6.0.1-MPI-Linux-Python3.12-x86_64/bin/pvbatch}" --force-offscreen-rendering "$@"
