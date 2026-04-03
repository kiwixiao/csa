#!/bin/bash
# =============================================================================
# Prepare subject folder for automated CSA pipeline
# =============================================================================
# Copies FFD motion data from MIRTK registration results to a clean
# subject folder structure. After running this, place your open-profile
# STL (exported from STAR-CCM+ or similar) in {subject}/surface/frame0.stl
#
# Usage:
#   ./prepare_subject.sh --reg-dir /path/to/registration_results --subject ENT001
#   ./prepare_subject.sh --reg-dir /path/to/ENT001_noalign_be0.0009 --subject ENT001 --output-dir /storage1
# =============================================================================

set -o errexit

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# --- CLI argument parsing ---
opt_reg_dir=""
opt_subject=""
opt_output_dir="."

while [[ $# -gt 0 ]]; do
    case $1 in
        --reg-dir)    opt_reg_dir="$2"; shift 2 ;;
        --subject)    opt_subject="$2"; shift 2 ;;
        --output-dir) opt_output_dir="$2"; shift 2 ;;
        --help)
            echo "Usage: prepare_subject.sh --reg-dir DIR --subject ID [--output-dir DIR]"
            echo ""
            echo "Options:"
            echo "  --reg-dir DIR      Path to MIRTK registration results (contains ffd_*.dof.gz)"
            echo "  --subject ID       Subject ID (e.g., ENT001)"
            echo "  --output-dir DIR   Where to create subject folder (default: current dir)"
            echo ""
            echo "After running, place your open-profile STL at:"
            echo "  {output-dir}/{subject}/surface/frame0.stl"
            exit 0
            ;;
        *) error "Unknown option: $1" ;;
    esac
done

# --- Validate inputs ---
[ -z "$opt_reg_dir" ] && error "Missing --reg-dir. Use --help for usage."
[ -z "$opt_subject" ] && error "Missing --subject. Use --help for usage."
[ ! -d "$opt_reg_dir" ] && error "Registration directory not found: $opt_reg_dir"

# Check for required files
ls "$opt_reg_dir"/ffd_*.dof.gz &>/dev/null || error "No ffd_*.dof.gz files found in $opt_reg_dir"

# --- Create subject folder structure ---
subject_dir="$opt_output_dir/$opt_subject"
info "Creating subject folder: $subject_dir"

mkdir -p "$subject_dir/registration"
mkdir -p "$subject_dir/surface"
mkdir -p "$subject_dir/branches"
mkdir -p "$subject_dir/motion/stl"
mkdir -p "$subject_dir/motion/centerlines"
mkdir -p "$subject_dir/csa"

# --- Copy FFD transforms ---
n_ffd=$(ls "$opt_reg_dir"/ffd_*.dof.gz 2>/dev/null | wc -l)
cp "$opt_reg_dir"/ffd_*.dof.gz "$subject_dir/registration/"
info "Copied $n_ffd FFD transforms"

# --- Copy ffds.csv ---
if [ -f "$opt_reg_dir/ffds.csv" ]; then
    cp "$opt_reg_dir/ffds.csv" "$subject_dir/registration/"
    info "Copied ffds.csv"
else
    warn "ffds.csv not found in $opt_reg_dir — you may need to create it manually"
fi

# --- Copy input.txt (timing parameters) ---
if [ -f "$opt_reg_dir/input.txt" ]; then
    cp "$opt_reg_dir/input.txt" "$subject_dir/registration/"
    info "Copied input.txt"
elif [ -f "$(dirname "$opt_reg_dir")/input.txt" ]; then
    cp "$(dirname "$opt_reg_dir")/input.txt" "$subject_dir/registration/"
    info "Copied input.txt (from parent dir)"
else
    warn "input.txt not found — you may need to create it manually"
fi

# --- Copy reference image (img_0 or seg_0) ---
ref_image=$(ls "$opt_reg_dir"/img_0.nii.gz "$opt_reg_dir"/staticImage_*_0.nii.gz "$opt_reg_dir"/seg_0.nii.gz 2>/dev/null | head -1)
if [ -n "$ref_image" ]; then
    cp "$ref_image" "$subject_dir/registration/img_0.nii.gz"
    info "Copied reference image: $(basename "$ref_image")"
else
    warn "No reference image (img_0.nii.gz) found"
fi

# --- Print summary ---
echo ""
echo "=========================================="
info "Subject folder ready: $subject_dir"
echo "=========================================="
echo ""
echo "Folder structure:"
echo "  $subject_dir/"
echo "    registration/        ← FFD motion data ($n_ffd transforms)"
echo "    surface/             ← PUT YOUR STL HERE"
echo "    branches/            ← branch_detector.py output (Step 2)"
echo "    motion/              ← interpolated time series (Step 3)"
echo "    csa/                 ← slicing results (Step 4)"
echo ""
echo -e "${YELLOW}NEXT STEP:${NC}"
echo "  1. Export your open-profile STL from STAR-CCM+ (with nostrils + trachea openings)"
echo "  2. Save it as: $subject_dir/surface/frame0.stl"
echo "  3. Run branch detection:"
echo "     conda run -n vmtk python python_slicer/branch_detector.py \\"
echo "       $subject_dir/surface/frame0.stl \\"
echo "       --subject-id $opt_subject \\"
echo "       --output-dir $subject_dir/branches"
echo ""
