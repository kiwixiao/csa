#!/usr/bin/env python3
"""
Advanced Airway Dynamics Analysis with Inter-Plane Volumes
===========================================================

Comprehensive analysis of upper airway dynamics during breathing cycle including:

1. Enhanced Shape Metrics (with explanations)
2. Centroid Movement Statistics (max, min, std per plane)
3. Airway Wall Movement Tracking (contour point tracking)
4. Inter-Plane Volume Analysis (volume between adjacent planes)
5. Volume Change Heatmaps (time vs band index)

Generates publication-ready PDF plots with embedded explanations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from pathlib import Path
import glob
import sys
import trimesh
from typing import List, Tuple, Dict
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

sys.path.insert(0, 'python_slicer')
from slicer.io_utils import read_vtk_centerline


def compute_circularity(area: float, perimeter: float) -> float:
    """
    Circularity = 4π × Area / Perimeter²
    1.0 = perfect circle, <1 = more elongated
    """
    if perimeter <= 0:
        return 0.0
    return (4.0 * np.pi * area) / (perimeter ** 2)


def compute_solidity(vertices_2d: np.ndarray, area: float) -> float:
    """
    Solidity = Area / ConvexHull_Area
    Measures concavity/irregularity (1.0 = convex, <1 = concave)
    """
    try:
        if len(vertices_2d) < 3:
            return 1.0
        hull = ConvexHull(vertices_2d)
        convex_area = hull.volume  # In 2D, 'volume' is actually area
        if convex_area <= 0:
            return 1.0
        return min(area / convex_area, 1.0)
    except Exception:
        return 1.0


def compute_eccentricity(major_axis: float, minor_axis: float) -> float:
    """
    Eccentricity = sqrt(1 - (b²/a²))
    where a = major axis, b = minor axis
    0 = circle, approaching 1 = very elongated
    """
    if major_axis <= 0 or minor_axis <= 0:
        return 0.0
    ratio = minor_axis / major_axis if major_axis > 0 else 1.0
    return np.sqrt(1 - ratio**2)


def load_plane_mesh(stl_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load plane mesh and return 3D vertices and 2D projection"""
    try:
        mesh = trimesh.load_mesh(stl_path)
        vertices_3d = mesh.vertices
        vertices_2d = vertices_3d[:, :2]  # Project to XY plane
        return vertices_3d, vertices_2d
    except Exception as e:
        print(f"  Warning: Could not load {stl_path}: {e}")
        return np.array([]), np.array([])


def compute_plane_normal(stl_path: str) -> np.ndarray:
    """Estimate plane normal from STL triangles"""
    try:
        mesh = trimesh.load_mesh(stl_path)
        normals = mesh.face_normals
        areas = mesh.area_faces
        if len(areas) == 0:
            return np.array([0, 0, 1])
        weighted_normal = np.average(normals, axis=0, weights=areas)
        norm = np.linalg.norm(weighted_normal)
        if norm > 0:
            return weighted_normal / norm
        return np.array([0, 0, 1])
    except Exception:
        return np.array([0, 0, 1])


def compute_interplane_volume(plane1_verts: np.ndarray, plane2_verts: np.ndarray) -> float:
    """
    Compute volume between two adjacent planes using convex hull approach.
    Creates a combined point cloud and computes convex hull volume.
    """
    try:
        if len(plane1_verts) < 3 or len(plane2_verts) < 3:
            return 0.0

        # Combine vertices from both planes
        combined_verts = np.vstack([plane1_verts, plane2_verts])

        # Compute convex hull in 3D
        hull = ConvexHull(combined_verts)

        return hull.volume
    except Exception as e:
        return 0.0


def track_contour_movement(plane_stl_files: List[str]) -> Dict:
    """
    Track individual contour point movement across time points.
    Uses nearest-neighbor matching between time points.

    Returns dict with:
      - max_point_displacement: maximum displacement of any point
      - mean_point_displacement: average displacement across all points
      - std_point_displacement: standard deviation
    """
    if len(plane_stl_files) < 2:
        return {"max_point_displacement": 0.0, "mean_point_displacement": 0.0, "std_point_displacement": 0.0}

    displacements = []

    for i in range(len(plane_stl_files) - 1):
        verts1, _ = load_plane_mesh(plane_stl_files[i])
        verts2, _ = load_plane_mesh(plane_stl_files[i + 1])

        if len(verts1) == 0 or len(verts2) == 0:
            continue

        # Find nearest neighbor matches
        distances = cdist(verts1, verts2, metric='euclidean')
        min_distances = distances.min(axis=1)

        displacements.extend(min_distances)

    if len(displacements) == 0:
        return {"max_point_displacement": 0.0, "mean_point_displacement": 0.0, "std_point_displacement": 0.0}

    displacements = np.array(displacements)

    return {
        "max_point_displacement": displacements.max(),
        "mean_point_displacement": displacements.mean(),
        "std_point_displacement": displacements.std()
    }


def analyze_enhanced_metrics(subject: str, partition: str, output_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze enhanced metrics from existing STL/CSV files.
    Returns: (enhanced_df, summary_df)
    """
    print("="*80)
    print("ADVANCED AIRWAY DYNAMICS ANALYSIS")
    print("="*80)
    print(f"Subject: {subject}")
    print(f"Partition: {partition}")
    print()

    # Load all Data CSV files
    csv_pattern = f"{partition}SlicedSTLs/*-Data.csv"
    csv_files = sorted(glob.glob(csv_pattern))

    if not csv_files:
        raise FileNotFoundError(f"No CSV files found matching: {csv_pattern}")

    print(f"Found {len(csv_files)} time points")

    # Load and combine all CSVs
    dfs = []
    for csv_file in csv_files:
        time_point = Path(csv_file).stem.replace('-Data', '')
        df = pd.read_csv(csv_file)
        df['time_point'] = time_point
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)

    print(f"Total rows: {len(combined_df)}")
    print(f"Unique planes: {combined_df['plane_index'].nunique()}")
    print()

    # Compute enhanced shape metrics
    print("Computing enhanced shape metrics...")
    enhanced_data = []

    for idx, row in combined_df.iterrows():
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx+1}/{len(combined_df)} rows...")

        # Circularity
        circularity = compute_circularity(row['area_mm2'], row['perimeter_mm'])

        # Eccentricity (from existing major/minor axes)
        eccentricity = compute_eccentricity(row['major_axis_mm'], row['minor_axis_mm'])

        # Load STL to compute solidity and normal
        stl_pattern = f"{partition}SlicedSTLs/{row['time_point']}-Planes-{row['plane_index']:03d}.stl"
        stl_files = glob.glob(stl_pattern)

        solidity = 1.0
        normal = np.array([0, 0, 1])

        if stl_files:
            _, vertices_2d = load_plane_mesh(stl_files[0])
            if len(vertices_2d) > 0:
                solidity = compute_solidity(vertices_2d, row['area_mm2'])
            normal = compute_plane_normal(stl_files[0])

        enhanced_data.append({
            'time_point': row['time_point'],
            'plane_index': row['plane_index'],
            'arc_length_mm': row['arc_length_mm'],
            'area_mm2': row['area_mm2'],
            'perimeter_mm': row['perimeter_mm'],
            'hydraulic_diameter_mm': row['hydraulic_diameter_mm'],
            'equivalent_diameter_mm': row['equivalent_diameter_mm'],
            'major_axis_mm': row['major_axis_mm'],
            'minor_axis_mm': row['minor_axis_mm'],
            'diameter_ratio': row['diameter_ratio'],
            'centroid_x': row['centroid_x'],
            'centroid_y': row['centroid_y'],
            'centroid_z': row['centroid_z'],
            'circularity': circularity,
            'eccentricity': eccentricity,
            'solidity': solidity,
            'normal_x': normal[0],
            'normal_y': normal[1],
            'normal_z': normal[2],
            'is_valid': row['is_valid']
        })

    enhanced_df = pd.DataFrame(enhanced_data)

    print("✓ Enhanced metrics computed")

    # Save enhanced metrics CSV
    output_csv = output_dir / f"{subject}_{partition}_enhanced_metrics.csv"
    enhanced_df.to_csv(output_csv, index=False)
    print(f"✓ Saved to: {output_csv.name}")
    print()

    # Compute summary statistics per plane
    print("Computing per-plane summary statistics...")
    summary_data = compute_summary_metrics(enhanced_df, partition)
    summary_df = pd.DataFrame(summary_data)

    # Save summary CSV
    summary_csv = output_dir / f"{subject}_{partition}_summary_metrics.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"✓ Saved summary to: {summary_csv.name}")
    print()

    return enhanced_df, summary_df


def compute_summary_metrics(enhanced_df: pd.DataFrame, partition: str) -> List[Dict]:
    """Compute summary metrics for each plane"""
    summary_data = []

    for plane_idx in sorted(enhanced_df['plane_index'].unique()):
        plane_data = enhanced_df[enhanced_df['plane_index'] == plane_idx]

        # CSA statistics
        area_mean = plane_data['area_mm2'].mean()
        area_std = plane_data['area_mm2'].std()
        area_min = plane_data['area_mm2'].min()
        area_max = plane_data['area_mm2'].max()
        area_range = area_max - area_min
        compliance = (area_range / area_mean) * 100 if area_mean > 0 else 0

        # Shape metrics
        circularity_mean = plane_data['circularity'].mean()
        circularity_std = plane_data['circularity'].std()
        eccentricity_mean = plane_data['eccentricity'].mean()
        solidity_mean = plane_data['solidity'].mean()

        # Centroid motion
        centroid_x = plane_data['centroid_x'].values
        centroid_y = plane_data['centroid_y'].values
        centroid_z = plane_data['centroid_z'].values

        # Displacements from mean position
        mean_x, mean_y, mean_z = centroid_x.mean(), centroid_y.mean(), centroid_z.mean()
        displacements_from_mean = np.sqrt(
            (centroid_x - mean_x)**2 + (centroid_y - mean_y)**2 + (centroid_z - mean_z)**2
        )

        # Sequential displacements (frame-to-frame)
        if len(centroid_x) > 1:
            displacements = np.sqrt(np.diff(centroid_x)**2 + np.diff(centroid_y)**2 + np.diff(centroid_z)**2)
            total_path_length = np.sum(displacements)
        else:
            displacements = np.array([0])
            total_path_length = 0.0

        # Rotation analysis (normal vector changes)
        normals = np.column_stack([
            plane_data['normal_x'].values,
            plane_data['normal_y'].values,
            plane_data['normal_z'].values
        ])

        max_rotation = 0.0
        for i in range(len(normals)):
            for j in range(i+1, len(normals)):
                dot_product = np.clip(np.dot(normals[i], normals[j]), -1.0, 1.0)
                angle = np.degrees(np.arccos(dot_product))
                max_rotation = max(max_rotation, angle)

        # Airway wall movement (contour points)
        stl_pattern = f"{partition}SlicedSTLs/*-Planes-{plane_idx:03d}.stl"
        stl_files = sorted(glob.glob(stl_pattern))
        wall_movement = track_contour_movement(stl_files)

        summary_data.append({
            'plane_index': plane_idx,
            'arc_length_mm': plane_data['arc_length_mm'].iloc[0],
            # CSA statistics
            'area_mean_mm2': area_mean,
            'area_std_mm2': area_std,
            'area_min_mm2': area_min,
            'area_max_mm2': area_max,
            'area_range_mm2': area_range,
            'area_cv_percent': (area_std / area_mean * 100) if area_mean > 0 else 0,
            'compliance_percent': compliance,
            # Shape metrics
            'circularity_mean': circularity_mean,
            'circularity_std': circularity_std,
            'eccentricity_mean': eccentricity_mean,
            'solidity_mean': solidity_mean,
            # 3D motion
            'total_path_length_mm': total_path_length,
            'max_displacement_mm': displacements_from_mean.max(),
            'mean_displacement_mm': displacements_from_mean.mean(),
            'std_displacement_mm': displacements_from_mean.std(),
            'max_rotation_deg': max_rotation,
            # Airway wall movement
            'wall_max_displacement_mm': wall_movement['max_point_displacement'],
            'wall_mean_displacement_mm': wall_movement['mean_point_displacement'],
            'wall_std_displacement_mm': wall_movement['std_point_displacement']
        })

    return summary_data


def compute_interplane_volumes(partition: str, enhanced_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute volumes between adjacent planes for all time points.
    Returns DataFrame with columns: band_index, time_point, volume_mm3
    """
    print("Computing inter-plane volumes...")

    time_points = sorted(enhanced_df['time_point'].unique())
    plane_indices = sorted(enhanced_df['plane_index'].unique())

    volume_data = []

    for t, time_point in enumerate(time_points):
        print(f"  [{t+1}/{len(time_points)}] {time_point}...", end=' ')

        volumes_computed = 0

        for i in range(len(plane_indices) - 1):
            plane1_idx = plane_indices[i]
            plane2_idx = plane_indices[i + 1]

            # Load both plane meshes
            stl1_pattern = f"{partition}SlicedSTLs/{time_point}-Planes-{plane1_idx:03d}.stl"
            stl2_pattern = f"{partition}SlicedSTLs/{time_point}-Planes-{plane2_idx:03d}.stl"

            stl1_files = glob.glob(stl1_pattern)
            stl2_files = glob.glob(stl2_pattern)

            if stl1_files and stl2_files:
                verts1, _ = load_plane_mesh(stl1_files[0])
                verts2, _ = load_plane_mesh(stl2_files[0])

                if len(verts1) > 0 and len(verts2) > 0:
                    volume = compute_interplane_volume(verts1, verts2)

                    volume_data.append({
                        'band_index': i,  # Band between plane i and i+1
                        'plane_pair': f"{plane1_idx}-{plane2_idx}",
                        'time_point': time_point,
                        'volume_mm3': volume
                    })

                    volumes_computed += 1

        print(f"{volumes_computed} bands")

    volume_df = pd.DataFrame(volume_data)
    print(f"✓ Computed {len(volume_df)} inter-plane volumes")
    print()

    return volume_df



def create_parameter_reference_page(pdf):
    """
    Create a comprehensive parameter definitions reference page.
    This serves as a 'legend' for all subsequent plots in the PDF.
    """
    fig = plt.figure(figsize=(11, 14))  # Letter size
    fig.suptitle('AIRWAY DYNAMICS ANALYSIS - PARAMETER REFERENCE',
                 fontsize=16, fontweight='bold', y=0.98)

    ax = fig.add_subplot(111)
    ax.axis('off')

    # Comprehensive parameter definitions text
    reference_text = """
GEOMETRIC PARAMETERS
────────────────────────────────────────────────────────────────────────────
• Cross-Sectional Area (CSA): Area of airway lumen at perpendicular plane (mm²)
  - Clinical significance: Directly correlates with airflow resistance
  - Interpretation: Smaller CSA → higher resistance, potential obstruction site

• Perimeter: Total boundary length of cross-sectional contour (mm)

• Plane Index: Position along airway centerline (sequential numbering)

• Band Index: Sequential numbering of airway segments 
  - Band N = volume between planes N and N+1

• Time Point (T0-TN): Breathing cycle phases (frame number in 4D imaging)
  - T0: Baseline/reference time point
  - Number of time points varies by subject/scan protocol


SHAPE METRICS
────────────────────────────────────────────────────────────────────────────
• Circularity: 4π × Area / Perimeter²
  - Range: 0 to 1.0 (1.0 = perfect circle, <1 = irregular)
  - Clinical: Drop in circularity = shape distortion (e.g., lateral collapse)

• Eccentricity: √(1 - (minor_axis/major_axis)²)
  - Range: 0 to 1 (0 = circle, →1 = elongated ellipse)
  - Clinical: Increase = elliptical deformation (common in soft palate region)

• Solidity: Contour_Area / Convex_Hull_Area
  - Range: 0 to 1.0 (1.0 = perfectly convex, <1 = irregular boundary)
  - Clinical: Decrease = boundary irregularities, may indicate turbulent flow


DYNAMIC METRICS (WALL MOVEMENT)
────────────────────────────────────────────────────────────────────────────
• Displacement: Euclidean distance (mm) between corresponding boundary points
  across time (point-to-point tracking)

• Max Displacement: Largest movement detected in cross-section (mm)
  - Clinical: High value = localized tissue mobility (potential collapse site)

• Mean Displacement: Average wall movement across all boundary points (mm)
  - Clinical: High value = overall airway compliance

• Std Displacement: Variability in wall movement across boundary points (mm)
  - Clinical: High std = non-uniform wall motion (complex deformation)

• Centroid Movement: Displacement of cross-section center of mass (mm)
  - Indicates airway shift/translation vs pure expansion/collapse


VOLUME METRICS
────────────────────────────────────────────────────────────────────────────
• Inter-Plane Volume: 3D airway lumen volume between adjacent cross-sections (mm³)
  - Formula: (Area_N + Area_N+1) / 2 × Distance_between_centroids
  - Clinical: Small volumes = narrow segments (potential obstruction sites)
            Large temporal variation = high compressibility/compliance
            Valleys across time = consistent narrowing (fixed restriction)

• Relative Volume Change: (Volume - T0_Volume) / T0_Volume × 100%
  - Purpose: Removes anatomical size bias by normalizing to baseline
  - Interpretation: 0% = no change; +% = expansion; -% = collapse
                  Large swings (±50%) = highly compliant/collapsible tissue


COMPLIANCE AND COLLAPSE RISK METRICS
────────────────────────────────────────────────────────────────────────────
• Coefficient of Variation (CV): Std(Volume) / Mean(Volume) × 100%
  - Measures temporal variability
  - High CV (>30%) = compliant tissue; Low CV (<10%) = stable rigid segment

• Min/Max Ratio: Min_Volume / Max_Volume
  - Range: 0 to 1 (0 = complete collapse, 1 = no change)
  - Low ratio (<0.5) = severe collapse potential

• Collapse Risk Score: CV × (1 - Min/Max_Ratio)
  - Combined metric indicating vulnerability to obstruction
  - High score = prime candidate for surgical/therapeutic intervention


VISUALIZATION ELEMENTS
────────────────────────────────────────────────────────────────────────────
• Heatmap Colors: Blue → Yellow (blue = smaller values, yellow = larger values)
  - Horizontal dark bands = consistently narrow regions (anatomical constriction)
  - Horizontal color variation = dynamic regions (collapse/expansion)
  - Vertical patterns = global airway response at specific breathing phases

• Box Plot Elements:
  - Box: Interquartile range (25th to 75th percentile)
  - Whiskers: Data range (excluding outliers)
  - Median Line: Red line inside box
  - Outliers: Individual points beyond whiskers
  - Narrow boxes = consistent values (stable anatomy)
  - Wide boxes = variable values (compliant/collapsible tissue)
"""

    ax.text(0.05, 0.95, reference_text, transform=ax.transAxes,
            fontsize=8.5, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    print("  ✓ Page 1: Parameter reference created")


def create_pdf_plots_with_explanations(enhanced_df: pd.DataFrame, summary_df: pd.DataFrame,
                                       volume_df: pd.DataFrame, output_dir: Path,
                                       subject: str, partition: str):
    """
    Create comprehensive PDF with all plots and embedded explanations.
    """
    print("Generating PDF plots with explanations...")

    pdf_path = output_dir / f"{subject}_{partition}_airway_dynamics_report.pdf"

    with PdfPages(pdf_path) as pdf:
        # Page 1: Parameter definitions reference
        create_parameter_reference_page(pdf)

        # Plot 1: CSA Dynamics
        create_csa_dynamics_plot_with_explanation(enhanced_df, pdf, partition)

        # Plot 2: Shape Metrics with explanation
        create_shape_metrics_plot_with_explanation(summary_df, pdf)

        # Plot 3: Centroid Movement (max, min, std bars)
        create_centroid_movement_bars(summary_df, pdf)

        # Plot 4: Airway Wall Movement
        create_wall_movement_plot(summary_df, pdf)

        # Plot 5: Inter-Plane Volume Dynamics (line plot like CSA)
        create_volume_dynamics_plot_with_explanation(volume_df, pdf)

        # Plot 6: Inter-Plane Volume Heatmap
        create_volume_heatmap(volume_df, pdf)

        # Plot 7: Relative Volume Changes (Normalized to T0)
        create_relative_volume_change_plot(volume_df, pdf)

        # Plot 9: Relative Volume Changes - Box Plot Distribution
        create_relative_volume_change_boxplot(volume_df, pdf)

        # Plot 10: Compliance/Collapse Risk
        create_compliance_plot_with_explanation(summary_df, pdf)

    print(f"✓ PDF report saved to: {pdf_path.name}")
    print()


def create_csa_dynamics_plot_with_explanation(enhanced_df, pdf, partition):
    """CSA dynamics over time - clean plot with proper title"""
    fig, ax1 = plt.subplots(figsize=(14, 10))

    # Get time points dynamically
    time_points = sorted(enhanced_df['time_point'].unique())
    num_timepoints = len(time_points)
    plane_indices = sorted(enhanced_df['plane_index'].unique())

    # Use colormap for time progression
    cmap = plt.cm.viridis
    colors = [cmap(i / len(time_points)) for i in range(len(time_points))]

    # Plot each time point as a separate line
    for i, time_point in enumerate(time_points):
        time_data = enhanced_df[enhanced_df['time_point'] == time_point]
        time_data_sorted = time_data.sort_values('plane_index')

        ax1.plot(time_data_sorted['plane_index'], time_data_sorted['area_mm2'],
                label=f'T{i}', alpha=0.6, linewidth=1.5, color=colors[i])

    ax1.set_xlabel('Plane Index (Along Airway)', fontsize=12)
    ax1.set_ylabel('Cross-Sectional Area (mm²)', fontsize=12)
    ax1.set_title(f'CSA Along Airway Length - {partition}', fontsize=14, fontweight='bold')
    ax1.legend(ncol=7, fontsize=8, title='Time Points')
    ax1.grid(True, alpha=0.3)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_shape_metrics_plot_with_explanation(summary_df, pdf):
    """Shape metrics along airway with explanation"""
    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(5, 1, height_ratios=[1.5, 3, 3, 3, 1.5], hspace=0.7)

    # Top explanation
    ax_text_top = fig.add_subplot(gs[0])
    ax_text_top.axis('off')
    explanation_top = (
        "SHAPE METRICS ALONG AIRWAY\n\n"
        "Parameter Definitions:\n"
        "• Circularity: 4π × Area / Perimeter² (1.0 = perfect circle, <1 = irregular)\n"
        "• Eccentricity: √(1 - (minor_axis/major_axis)²) (0 = circle, →1 = elongated ellipse)\n"
        "• Solidity: Contour_Area / Convex_Hull_Area (1.0 = perfectly convex, <1 = irregular boundary)\n"
        "• Plane Index: Position along airway centerline\n\n"
        "These metrics quantify the geometry of each cross-section:"
    )
    ax_text_top.text(0.05, 0.5, explanation_top, transform=ax_text_top.transAxes,
                    fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Circularity
    ax1 = fig.add_subplot(gs[1])
    ax1.plot(summary_df['plane_index'], summary_df['circularity_mean'], 'o-', color='blue', linewidth=2)
    ax1.fill_between(summary_df['plane_index'],
                     summary_df['circularity_mean'] - summary_df['circularity_std'],
                     summary_df['circularity_mean'] + summary_df['circularity_std'],
                     alpha=0.3, color='blue')
    ax1.set_ylabel('Circularity', fontsize=11)
    ax1.set_title('Circularity: 1.0 = perfect circle, <1 = elongated', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1.05])

    # Eccentricity
    ax2 = fig.add_subplot(gs[2])
    ax2.plot(summary_df['plane_index'], summary_df['eccentricity_mean'], 'o-', color='green', linewidth=2)
    ax2.set_ylabel('Eccentricity', fontsize=11)
    ax2.set_title('Eccentricity: 0 = circle, approaching 1 = very elongated', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1.05])

    # Solidity
    ax3 = fig.add_subplot(gs[3])
    ax3.plot(summary_df['plane_index'], summary_df['solidity_mean'], 'o-', color='orange', linewidth=2)
    ax3.set_xlabel('Plane Index (Along Airway)', fontsize=11)
    ax3.set_ylabel('Solidity', fontsize=11)
    ax3.set_title('Solidity: 1.0 = convex, <1 = concave/irregular', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim([0, 1.05])

    # Bottom explanation
    ax_text_bottom = fig.add_subplot(gs[4])
    ax_text_bottom.axis('off')
    explanation_bottom = (
        "Clinical Significance:\n"
        "• Circularity drop = shape distortion (e.g., lateral collapse)\n"
        "• Eccentricity increase = elliptical deformation\n"
        "• Solidity decrease = boundary irregularities or indentations\n"
        "• Low circularity/high eccentricity → elliptical cross-section (common in soft palate region)\n"
        "• Low solidity → irregular/concave shape (may indicate turbulent flow regions)"
    )
    ax_text_bottom.text(0.05, 0.5, explanation_bottom, transform=ax_text_bottom.transAxes,
                       fontsize=9, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_centroid_movement_bars(summary_df, pdf):
    """Bar chart showing max, min, std of centroid displacement for each plane"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    x = summary_df['plane_index']
    width = 0.6

    # Max displacement
    ax1.bar(x, summary_df['max_displacement_mm'], width=width, color='darkred', alpha=0.7, label='Max')
    ax1.set_ylabel('Max Displacement (mm)', fontsize=11)
    ax1.set_title('Maximum 3D Centroid Displacement from Mean Position', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend()

    # Mean ± Std displacement
    ax2.bar(x, summary_df['mean_displacement_mm'], width=width, color='steelblue', alpha=0.7,
            yerr=summary_df['std_displacement_mm'], capsize=3, label='Mean ± Std')
    ax2.set_xlabel('Plane Index (Along Airway)', fontsize=11)
    ax2.set_ylabel('Mean Displacement (mm)', fontsize=11)
    ax2.set_title('Mean 3D Centroid Displacement with Standard Deviation', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend()

    fig.text(0.5, 0.02, 'Higher values indicate greater tissue mobility during breathing',
            ha='center', fontsize=10, style='italic')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_wall_movement_plot(summary_df, pdf):
    """Airway wall movement (contour point tracking)"""
    fig = plt.figure(figsize=(14, 14))
    gs = GridSpec(5, 1, height_ratios=[1.5, 3, 3, 3, 1.5], hspace=0.6)

    # Top explanation
    ax_text_top = fig.add_subplot(gs[0])
    ax_text_top.axis('off')
    explanation_top = (
        "AIRWAY WALL MOVEMENT ANALYSIS\n\n"
        "Parameter Definitions:\n"
        "• Displacement: Euclidean distance (mm) between corresponding boundary points across time\n"
        "• Max Displacement: Largest movement detected in the cross-section (mm)\n"
        "• Mean Displacement: Average wall movement across all boundary points (mm)\n"
        "• Std Displacement: Variability in wall movement (mm)\n"
        "• Plane Index: Position along airway centerline"
    )
    ax_text_top.text(0.05, 0.5, explanation_top, transform=ax_text_top.transAxes,
                    fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Plot axes
    axes = [fig.add_subplot(gs[1]), fig.add_subplot(gs[2]), fig.add_subplot(gs[3])]

    x = summary_df['plane_index']
    width = 0.6

    # Max wall displacement
    axes[0].bar(x, summary_df['wall_max_displacement_mm'], width=width, color='darkgreen', alpha=0.7)
    axes[0].set_ylabel('Max Wall Displacement (mm)', fontsize=11)
    axes[0].set_title('Maximum Airway Wall Point Displacement', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y')

    # Mean wall displacement
    axes[1].bar(x, summary_df['wall_mean_displacement_mm'], width=width, color='mediumseagreen', alpha=0.7)
    axes[1].set_ylabel('Mean Wall Displacement (mm)', fontsize=11)
    axes[1].set_title('Mean Airway Wall Point Displacement', fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y')

    # Std wall displacement
    axes[2].bar(x, summary_df['wall_std_displacement_mm'], width=width, color='lightseagreen', alpha=0.7)
    axes[2].set_xlabel('Plane Index (Along Airway)', fontsize=11)
    axes[2].set_ylabel('Std Wall Displacement (mm)', fontsize=11)
    axes[2].set_title('Standard Deviation of Wall Point Displacement', fontsize=12, fontweight='bold')
    axes[2].grid(True, alpha=0.3, axis='y')

    # Bottom explanation
    ax_text_bottom = fig.add_subplot(gs[4])
    ax_text_bottom.axis('off')
    explanation_bottom = (
        "Clinical Interpretation:\n"
        "• High max displacement = localized tissue mobility (potential collapse site)\n"
        "• High mean displacement = overall airway compliance\n"
        "• High std = non-uniform wall motion (complex deformation patterns)"
    )
    ax_text_bottom.text(0.05, 0.5, explanation_bottom, transform=ax_text_bottom.transAxes,
                       fontsize=9, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_volume_dynamics_plot_with_explanation(volume_df, pdf):
    """Volume dynamics along airway with embedded explanation (same format as CSA)"""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 1, height_ratios=[1.5, 3, 3, 3], hspace=0.6)

    # Explanation text box
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis('off')
    explanation = (
        "INTER-PLANE VOLUME DYNAMICS\n\n"
        "What it shows: Spatial distribution of airway segment volumes along the airway length at different breathing phases.\n\n"
        "Parameter Definitions:\n"
        "• Band Index: Sequential numbering of airway segments (band N = volume between planes N and N+1)\n"
        "• Inter-Plane Volume: 3D airway lumen volume between adjacent cross-sections (mm³)\n"
        "  Computed as: (Area_N + Area_N+1) / 2 × Distance_between_centroids\n"
        f"• Time Point: Breathing cycle phases ({len(sorted(volume_df['time_point'].unique()))} total frames)\n\n"
        "Interpretation:\n"
        "• X-axis = Band index (airway segment position between adjacent planes)\n"
        "• Each line = One time point (breathing phase)\n"
        "• Large variation between lines = high compliance/collapsibility\n"
        "• Small volumes (dips) = potential narrowing/collapse sites"
    )
    ax_text.text(0.05, 0.5, explanation, transform=ax_text.transAxes,
                fontsize=10, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.5))

    # Plot data in remaining space
    ax1 = fig.add_subplot(gs[1:])

    time_points = sorted(volume_df['time_point'].unique())
    band_indices = sorted(volume_df['band_index'].unique())

    # Use colormap for time progression
    cmap = plt.cm.viridis
    colors = [cmap(i / len(time_points)) for i in range(len(time_points))]

    # Plot each time point as a separate line
    for i, time_point in enumerate(time_points):
        time_data = volume_df[volume_df['time_point'] == time_point]
        time_data_sorted = time_data.sort_values('band_index')

        ax1.plot(time_data_sorted['band_index'], time_data_sorted['volume_mm3'],
                label=f'T{i}', alpha=0.6, linewidth=1.5, color=colors[i])

    ax1.set_xlabel('Band Index (Along Airway)', fontsize=12)
    ax1.set_ylabel('Volume (mm³)', fontsize=12)
    ax1.set_title('Inter-Plane Volume Along Airway Length', fontsize=14, fontweight='bold')
    ax1.legend(ncol=7, fontsize=8, title='Time Points')
    ax1.grid(True, alpha=0.3)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_temporal_volume_comparison(volume_df, pdf, bands_to_plot=[69, 70, 71, 72]):
    """
    Show how specific bands' volumes change throughout the breathing cycle.
    This helps distinguish between consistently wide bands vs breathing-induced changes.
    """
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 1, height_ratios=[1.5, 3, 3, 3], hspace=0.6)

    # Explanation text box
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis('off')
    explanation = (
        "TEMPORAL VOLUME COMPARISON - SELECTED BANDS\n\n"
        "What it shows: How specific airway segment volumes change throughout the breathing cycle.\n\n"
        "Parameter Definitions:\n"
        "• Time Point: Breathing cycle phase (0-20, typically representing full respiratory cycle)\n"
        "• Volume (mm³): Inter-plane volume for each band\n"
        "• Band Index: Specific airway segment (between planes N and N+1)\n\n"
        "Interpretation:\n"
        "• Peaks = expansion phases (inspiration or muscle activation)\n"
        "• Valleys = collapse phases (expiration or muscle relaxation)\n"
        "• Peak-to-valley amplitude = degree of compliance/collapsibility\n"
        "• Flat lines = rigid segments (bone/cartilage-supported)\n"
        "• Asymmetric patterns = complex breathing dynamics"
    )
    ax_text.text(0.05, 0.5, explanation, transform=ax_text.transAxes,
                fontsize=10, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.5))

    # Plot data in remaining space
    ax1 = fig.add_subplot(gs[1:])

    time_points = sorted(volume_df['time_point'].unique())

    # Filter for bands that exist in the data
    available_bands = volume_df['band_index'].unique()
    bands_to_plot = [b for b in bands_to_plot if b in available_bands]

    if len(bands_to_plot) == 0:
        ax1.text(0.5, 0.5, 'No data for requested bands',
                ha='center', va='center', fontsize=14)
    else:
        # Use distinct colors for each band
        colors = plt.cm.tab10(np.linspace(0, 1, len(bands_to_plot)))

        for i, band_idx in enumerate(bands_to_plot):
            band_data = volume_df[volume_df['band_index'] == band_idx].copy()
            band_data = band_data.sort_values('time_point')

            # Extract time point numbers for X-axis
            time_indices = [time_points.index(tp) for tp in band_data['time_point']]
            volumes = band_data['volume_mm3'].values

            # Plot with markers and lines
            ax1.plot(time_indices, volumes,
                    marker='o', markersize=6, linewidth=2,
                    label=f'Band {band_idx}', color=colors[i], alpha=0.8)

            # Add mean line for reference
            mean_vol = volumes.mean()
            ax1.axhline(mean_vol, color=colors[i], linestyle='--',
                       linewidth=1, alpha=0.3)

        ax1.set_xlabel('Time Point Index (Breathing Phase)', fontsize=12)
        ax1.set_ylabel('Volume (mm³)', fontsize=12)
        ax1.set_title('Temporal Volume Evolution - Selected Bands', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10, loc='best')
        ax1.grid(True, alpha=0.3)

        # Set integer X-axis ticks
        ax1.set_xticks(range(len(time_points)))
        ax1.set_xticklabels([f'T{i}' for i in range(len(time_points))], rotation=45)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_relative_volume_change_plot(volume_df, pdf):
    """
    Show relative volume changes (% change from frame 0) across all bands.
    X-axis = band index, Y-axis = volume %, each line = different time point (like CSA plot).
    This removes the bias from bands with different absolute widths.
    """
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 1, height_ratios=[1.5, 3, 3, 3], hspace=0.6)

    # Explanation text box
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis('off')
    explanation = (
        "RELATIVE VOLUME CHANGES (NORMALIZED TO FRAME 0)\n\n"
        "What it shows: Percentage change in volume relative to the first time point (T0).\n\n"
        "Parameter Definitions:\n"
        "• Relative Volume Change: (Volume - T0_Volume) / T0_Volume × 100%\n"
        "  This removes anatomical size bias by normalizing to baseline\n"
        "• Band Index: Position along airway (each band = volume between 2 planes)\n"
        f"• Time Point (T0-T{len(sorted(volume_df['time_point'].unique()))-1}): Breathing cycle phases ({len(sorted(volume_df['time_point'].unique()))} total frames)\n"
        "• T0 = Baseline/reference time point\n\n"
        "Interpretation:\n"
        "• 0% = no change from baseline\n"
        "• Positive values (+%) = expansion from baseline\n"
        "• Negative values (-%) = collapse/compression from baseline\n"
        "• Large swings (e.g., ±50%) = highly compliant/collapsible tissue"
    )
    ax_text.text(0.05, 0.5, explanation, transform=ax_text.transAxes,
                fontsize=10, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

    # Plot data in remaining space
    ax1 = fig.add_subplot(gs[1:])

    # Get all time points and bands
    time_points = sorted(volume_df['time_point'].unique())
    all_bands = sorted(volume_df['band_index'].unique())

    if len(time_points) == 0 or len(all_bands) == 0:
        ax1.text(0.5, 0.5, 'No volume data available',
                ha='center', va='center', fontsize=14)
    else:
        # Use colormap for time points (gradient showing progression through breathing cycle)
        colors = plt.cm.viridis(np.linspace(0, 1, len(time_points)))

        # Loop over each time point
        for i, time_point in enumerate(time_points):
            # Get all bands for this time point
            time_data = volume_df[volume_df['time_point'] == time_point].copy()
            time_data = time_data.sort_values('band_index')

            band_indices = time_data['band_index'].values
            volumes = time_data['volume_mm3'].values

            # Compute relative changes compared to T0
            # For each band, find its baseline volume at T0
            relative_changes = []
            valid_band_indices = []

            for band_idx, vol in zip(band_indices, volumes):
                # Get T0 volume for this band
                t0_data = volume_df[(volume_df['band_index'] == band_idx) &
                                   (volume_df['time_point'] == time_points[0])]

                if not t0_data.empty:
                    baseline_volume = t0_data.iloc[0]['volume_mm3']
                    if baseline_volume > 0:
                        rel_change = ((vol - baseline_volume) / baseline_volume) * 100.0
                        relative_changes.append(rel_change)
                        valid_band_indices.append(band_idx)

            # Plot this time point
            if len(relative_changes) > 0:
                time_idx = time_points.index(time_point)
                label = f'T{time_idx}' if i % 3 == 0 or i == len(time_points) - 1 else None  # Label every 3rd to avoid clutter

                ax1.plot(valid_band_indices, relative_changes,
                        linewidth=1.5 if i == 0 else 1.0,  # Make T0 (baseline) thicker
                        label=label, color=colors[i], alpha=0.7)

        # Add zero reference line
        ax1.axhline(0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5, label='Baseline (T0)')

        ax1.set_xlabel('Band Index (Position Along Airway)', fontsize=12)
        ax1.set_ylabel('Volume Change from T0 (%)', fontsize=12)
        ax1.set_title('Relative Volume Changes Along Airway Length Throughout Breathing Cycle',
                     fontsize=14, fontweight='bold')
        ax1.legend(fontsize=9, loc='best', ncol=2)
        ax1.grid(True, alpha=0.3)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_relative_volume_change_boxplot(volume_df, pdf):
    """
    Box plot showing distribution of relative volume changes for each band.
    X-axis = band index, Y-axis = volume %, each box = distribution across all time points.
    """
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 1, height_ratios=[1.5, 3, 3, 1.5], hspace=0.6)

    # Explanation text box
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis('off')
    explanation = (
        "RELATIVE VOLUME CHANGE DISTRIBUTION (BOX PLOT)\n\n"
        "What it shows: Statistical distribution of volume changes for each band across all time points.\n\n"
        "Parameter Definitions:\n"
        "• Band Index: Sequential numbering of airway segments along length\n"
        "• Relative Volume Change: (Volume - T0_Volume) / T0_Volume × 100%\n"
        "• Box: Shows interquartile range (25th to 75th percentile)\n"
        "• Whiskers: Show data range (excluding outliers)\n"
        "• Median Line: Red line inside box showing median value\n"
        "• Outliers: Individual points beyond whiskers"
    )
    ax_text.text(0.05, 0.5, explanation, transform=ax_text.transAxes,
                fontsize=10, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

    # Plot data in remaining space
    ax1 = fig.add_subplot(gs[1:3])

    # Get all time points and bands
    time_points = sorted(volume_df['time_point'].unique())
    all_bands = sorted(volume_df['band_index'].unique())

    if len(time_points) == 0 or len(all_bands) == 0:
        ax1.text(0.5, 0.5, 'No volume data available',
                ha='center', va='center', fontsize=14)
    else:
        # Prepare data for box plot
        # For each band, collect all relative volume changes across time points
        box_data = []
        box_positions = []

        for band_idx in all_bands:
            # Get all time points for this band
            band_data = volume_df[volume_df['band_index'] == band_idx].copy()
            band_data = band_data.sort_values('time_point')

            # Get T0 volume
            t0_data = band_data[band_data['time_point'] == time_points[0]]
            if not t0_data.empty:
                baseline_volume = t0_data.iloc[0]['volume_mm3']

                if baseline_volume > 0:
                    # Compute relative changes for all time points
                    volumes = band_data['volume_mm3'].values
                    relative_changes = ((volumes - baseline_volume) / baseline_volume) * 100.0

                    box_data.append(relative_changes)
                    box_positions.append(band_idx)

        # Create box plot
        if len(box_data) > 0:
            bp = ax1.boxplot(box_data, positions=box_positions, widths=0.6,
                            patch_artist=True, showfliers=True,
                            medianprops=dict(color='red', linewidth=2),
                            boxprops=dict(facecolor='lightblue', alpha=0.7),
                            whiskerprops=dict(linewidth=1.5),
                            capprops=dict(linewidth=1.5))

            # Add zero reference line
            ax1.axhline(0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

            ax1.set_xlabel('Band Index (Position Along Airway)', fontsize=12)
            ax1.set_ylabel('Volume Change from T0 (%)', fontsize=12)
            ax1.set_title('Distribution of Relative Volume Changes Across All Time Points',
                         fontsize=14, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y')

    # Bottom explanation
    ax_text_bottom = fig.add_subplot(gs[3])
    ax_text_bottom.axis('off')
    explanation_bottom = (
        "Clinical Interpretation:\n"
        "• Narrow boxes = consistent volume (stable anatomy)\n"
        "• Wide boxes = variable volume (compliant/collapsible tissue)\n"
        "• Box position above 0 = predominantly expansion; below 0 = predominantly collapse"
    )
    ax_text_bottom.text(0.05, 0.5, explanation_bottom, transform=ax_text_bottom.transAxes,
                       fontsize=9, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_volume_heatmap(volume_df, pdf):
    """Heatmap of inter-plane volumes: band index vs time"""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(4, 1, height_ratios=[1, 3, 0.3, 1], hspace=0.6)

    # Prepare heatmap data (get time points first for explanation text)
    band_indices = sorted(volume_df['band_index'].unique())
    time_points = sorted(volume_df['time_point'].unique())
    num_timepoints = len(time_points)

    # Explanation
    ax_text = fig.add_subplot(gs[0])
    ax_text.axis('off')
    explanation = (
        "INTER-PLANE VOLUME DYNAMICS\n\n"
        "What it shows: Volume of airway segments between adjacent planes throughout breathing cycle.\n\n"
        "Parameter Definitions:\n"
        "• Band Index (X-axis): Position along airway (each band = volume between 2 planes)\n"
        f"• Time Point (Y-axis): Breathing cycle progression (T0-T{num_timepoints-1}, {num_timepoints} frames total)\n"
        "• Color Intensity: Inter-plane volume magnitude (mm³)\n"
        "  Blue = smaller volumes, Yellow = larger volumes\n\n"
        "Interpretation:\n"
        "• Horizontal dark bands = consistently narrow regions (anatomical constriction)\n"
        "• Horizontal color variation = dynamic regions (collapse/expansion)\n"
        "• Vertical patterns = global airway response at specific breathing phases\n"
        "• Localized dark spots = transient collapse events\n"
        "• Vertical patterns (same band, different times) = breathing-induced volume changes\n"
        "• Horizontal patterns = anatomical regions with similar compliance"
    )
    ax_text.text(0.05, 0.5, explanation, transform=ax_text.transAxes,
                fontsize=10, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.5))

    # Create matrix: rows = time points, columns = band indices
    volume_matrix = np.zeros((len(time_points), len(band_indices)))

    for i, time_point in enumerate(time_points):
        for j, band_idx in enumerate(band_indices):
            vol_data = volume_df[(volume_df['band_index'] == band_idx) &
                                (volume_df['time_point'] == time_point)]
            if len(vol_data) > 0:
                volume_matrix[i, j] = vol_data.iloc[0]['volume_mm3']

    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1])
    im = ax_heatmap.imshow(volume_matrix, aspect='auto', cmap='YlOrRd', interpolation='nearest')

    ax_heatmap.set_xlabel('Band Index (Along Airway)', fontsize=12)
    ax_heatmap.set_ylabel('Time Point (Breathing Cycle)', fontsize=12)
    ax_heatmap.set_title('Inter-Plane Volume Distribution Along Airway', fontsize=14, fontweight='bold')

    # Colorbar
    ax_cbar = fig.add_subplot(gs[2])
    cbar = plt.colorbar(im, cax=ax_cbar, orientation='horizontal')
    cbar.set_label('Volume (mm³)', fontsize=11)

    # Bottom explanation
    ax_text_bottom = fig.add_subplot(gs[3])
    ax_text_bottom.axis('off')
    explanation_bottom = (
        "Clinical Significance:\n"
        "• Vertical color variation (same band, different times) → compliant/collapsible region\n"
        "• Consistent colors vertically → rigid region (bony structures)\n"
        "• Overall pattern indicates breathing mechanics and collapse risk at each anatomical location"
    )
    ax_text_bottom.text(0.05, 0.5, explanation_bottom, transform=ax_text_bottom.transAxes,
                       fontsize=9, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def create_compliance_plot_with_explanation(summary_df, pdf):
    """Compliance/collapse risk with explanation"""
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 1, height_ratios=[1, 3, 1], hspace=0.6)

    # Top explanation
    ax_text_top = fig.add_subplot(gs[0])
    ax_text_top.axis('off')
    explanation_top = (
        "COMPLIANCE AND COLLAPSE RISK\n\n"
        "Parameter Definitions:\n"
        "• Compliance: (Max_CSA - Min_CSA) / Mean_CSA × 100%\n"
        "  Measures temporal variability in airway area (higher = more dynamic/compliant)\n"
        "• Plane Index: Position along airway centerline\n"
        "• Max CSA, Min CSA: Maximum and minimum cross-sectional areas across all time points (mm²)\n"
        "• Mean CSA: Average cross-sectional area across all time points (mm²)\n\n"
        "Interpretation:\n"
        "• >30% (red) = High collapse risk, very compliant tissue\n"
        "• 15-30% (orange) = Moderate compliance\n"
        "• <15% (green) = Low compliance, stable structure"
    )
    ax_text_top.text(0.05, 0.5, explanation_top, transform=ax_text_top.transAxes,
                    fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Compliance plot
    ax = fig.add_subplot(gs[1])

    colors = ['green' if c < 15 else 'orange' if c < 30 else 'red'
             for c in summary_df['compliance_percent']]

    ax.bar(summary_df['plane_index'], summary_df['compliance_percent'],
          color=colors, alpha=0.7, width=0.6)

    ax.axhline(y=15, color='orange', linestyle='--', linewidth=2, label='15% threshold')
    ax.axhline(y=30, color='red', linestyle='--', linewidth=2, label='30% threshold')

    ax.set_xlabel('Plane Index (Along Airway)', fontsize=12)
    ax.set_ylabel('Compliance (%)', fontsize=12)
    ax.set_title('Airway Compliance and Collapse Risk', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    # Bottom explanation
    ax_text_bottom = fig.add_subplot(gs[2])
    ax_text_bottom.axis('off')
    explanation_bottom = (
        "Clinical Significance:\n"
        "High compliance regions are most likely to collapse during sleep, contributing to obstructive sleep apnea (OSA).\n"
        "Regions consistently above 30% compliance may benefit from surgical intervention or positive airway pressure therapy."
    )
    ax_text_bottom.text(0.05, 0.5, explanation_bottom, transform=ax_text_bottom.transAxes,
                       fontsize=9, verticalalignment='center',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Advanced Airway Dynamics Analysis')
    parser.add_argument('subject', help='Subject ID (e.g., OSAMRI037)')
    parser.add_argument('partition', help='Partition name (e.g., LeftNoseDecending)')
    args = parser.parse_args()

    # Output directory
    output_dir = Path('.')

    # Run analysis
    enhanced_df, summary_df = analyze_enhanced_metrics(args.subject, args.partition, output_dir)

    # Compute inter-plane volumes
    print("Computing inter-plane volumes (this may take a few minutes)...")
    volume_df = compute_interplane_volumes(args.partition, enhanced_df)

    # Save volume data
    volume_csv = output_dir / f"{args.subject}_{args.partition}_interplane_volumes.csv"
    volume_df.to_csv(volume_csv, index=False)
    print(f"✓ Saved inter-plane volumes to: {volume_csv.name}")
    print()

    # Create PDF plots
    create_pdf_plots_with_explanations(enhanced_df, summary_df, volume_df,
                                       output_dir, args.subject, args.partition)

    print("="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"Enhanced metrics: {args.subject}_{args.partition}_enhanced_metrics.csv")
    print(f"Summary metrics: {args.subject}_{args.partition}_summary_metrics.csv")
    print(f"Inter-plane volumes: {args.subject}_{args.partition}_interplane_volumes.csv")
    print(f"PDF report: {args.subject}_{args.partition}_airway_dynamics_report.pdf")
    print("="*80)


if __name__ == "__main__":
    main()
