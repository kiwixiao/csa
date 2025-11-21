#!/bin/bash
# Setup conda environment for Python Slicer

set -e  # Exit on error

echo "======================================================================"
echo "Python Slicer - Conda Environment Setup"
echo "======================================================================"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install Anaconda or Miniconda first."
    exit 1
fi

echo ""
echo "Creating conda environment 'csa'..."
conda env create -f environment.yml

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================================"
    echo "Setup complete!"
    echo "======================================================================"
    echo ""
    echo "To activate the environment:"
    echo "  conda activate csa"
    echo ""
    echo "To run the slicer from project root:"
    echo "  cd /storage1/Qiwei/CSATemplate"
    echo "  conda activate csa"
    echo "  python python_slicer/main.py OSAMRI037 LeftNoseDecending"
    echo ""
    echo "Or use the convenience script:"
    echo "  ./run_slicer.sh OSAMRI037 LeftNoseDecending"
    echo ""
else
    echo ""
    echo "Environment creation failed. If 'csa' already exists, remove it first:"
    echo "  conda env remove -n csa"
    echo "Then run this script again."
    exit 1
fi
