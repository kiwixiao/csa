#!/bin/bash
# Setup script for Python Slicer

set -e  # Exit on error

echo "======================================================================"
echo "Python Slicer - Setup Script"
echo "======================================================================"

# Check Python version
echo ""
echo "Checking Python version..."
python3 --version

if [ $? -ne 0 ]; then
    echo "Error: Python 3 not found. Please install Python 3.8 or later."
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo ""
    echo "Creating virtual environment..."
    python3 -m venv venv
else
    echo ""
    echo "Virtual environment already exists."
fi

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo ""
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo ""
echo "Installing dependencies..."
pip install -r requirements.txt

# Verify installation
echo ""
echo "Verifying installation..."
python -c "
import numpy
import scipy
import trimesh
import pyvista
import shapely
import pandas
print('✓ All core dependencies installed successfully!')
"

if [ $? -eq 0 ]; then
    echo ""
    echo "======================================================================"
    echo "Setup complete!"
    echo "======================================================================"
    echo ""
    echo "To activate the virtual environment in the future, run:"
    echo "  source venv/bin/activate"
    echo ""
    echo "To run the slicer:"
    echo "  python main.py <subject_name> <partition>"
    echo ""
    echo "Example:"
    echo "  python main.py DYMOSA801 LeftNoseDecending"
    echo ""
else
    echo ""
    echo "Error: Some dependencies failed to install."
    echo "Please check the error messages above."
    exit 1
fi
