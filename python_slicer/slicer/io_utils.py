"""
I/O Utilities for STL and VTK file operations
Replaces: STL_Import.m, read_vtk.m, csvwrite_with_headers.m
"""

import numpy as np
import trimesh
import pyvista as pv
import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, Any


def read_stl(filepath: str) -> trimesh.Trimesh:
    """
    Read STL file (ASCII or binary)

    Replaces: STL_Import.m

    Args:
        filepath: Path to STL file

    Returns:
        trimesh.Trimesh object with vertices, faces, and normals
    """
    try:
        mesh = trimesh.load_mesh(filepath)

        # Ensure mesh is valid
        if not isinstance(mesh, trimesh.Trimesh):
            raise ValueError(f"Loaded mesh is not a valid Trimesh object: {type(mesh)}")

        # Check for empty mesh
        if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
            raise ValueError(f"Empty mesh loaded from {filepath}")

        print(f"Loaded STL: {Path(filepath).name}")
        print(f"  Vertices: {len(mesh.vertices)}, Faces: {len(mesh.faces)}")

        return mesh

    except Exception as e:
        raise IOError(f"Failed to read STL file {filepath}: {str(e)}")


def read_vtk_centerline(filepath: str) -> np.ndarray:
    """
    Read VTK file containing centerline points

    Replaces: read_vtk.m

    Args:
        filepath: Path to VTK file

    Returns:
        np.ndarray of shape (N, 3) with centerline points
    """
    try:
        # Read VTK file using pyvista
        mesh = pv.read(filepath)

        # Extract points (vertices)
        points = np.array(mesh.points)

        # Get unique points maintaining order (replaces unique(pos','rows','stable'))
        # This is important to maintain the centerline order
        _, unique_indices = np.unique(points, axis=0, return_index=True)
        unique_indices = np.sort(unique_indices)  # Maintain original order
        points = points[unique_indices]

        if len(points) < 2:
            raise ValueError(f"Centerline has fewer than 2 points: {len(points)}")

        print(f"Loaded VTK centerline: {Path(filepath).name}")
        print(f"  Points: {len(points)}")

        return points

    except Exception as e:
        raise IOError(f"Failed to read VTK file {filepath}: {str(e)}")


def write_cross_section_stl(filepath: str, vertices: np.ndarray, faces: np.ndarray) -> None:
    """
    Write cross-section mesh to STL file

    Replaces: stlwrite.m

    Args:
        filepath: Output STL file path
        vertices: Nx3 array of vertex coordinates
        faces: Mx3 array of triangle indices
    """
    try:
        # Create trimesh object
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

        # Ensure output directory exists
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        # Export to STL
        mesh.export(filepath)

    except Exception as e:
        raise IOError(f"Failed to write STL file {filepath}: {str(e)}")


def write_measurements_csv(filepath: str, data: Dict[str, np.ndarray],
                           flip_order: bool = True) -> None:
    """
    Write measurement data to CSV file

    Replaces: csvwrite and custom CSV writing in MATLAB

    Args:
        filepath: Output CSV file path
        data: Dictionary with column names as keys and arrays as values
              Expected keys: 'arc_length', 'area', 'hydraulic_diameter', etc.
        flip_order: If True, reverse the order of rows (MATLAB compatibility)
    """
    try:
        # Create DataFrame
        df = pd.DataFrame(data)

        # Flip order if requested (MATLAB code does dataToWrite(:,1) = arcLength(end:-1:1)')
        if flip_order:
            df = df.iloc[::-1].reset_index(drop=True)

        # Ensure output directory exists
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        # Write to CSV
        df.to_csv(filepath, index=False)

        print(f"Wrote CSV: {Path(filepath).name} ({len(df)} rows)")

    except Exception as e:
        raise IOError(f"Failed to write CSV file {filepath}: {str(e)}")


def write_all_planes_stl(filepath: str, all_vertices: list, all_faces: list) -> None:
    """
    Write all cross-section planes combined into a single STL file

    Replaces: stlwrite(['./' STLfilename(1:end-4) '-Planes-All' '.stl'], tri, nodes);

    Args:
        filepath: Output STL file path
        all_vertices: List of vertex arrays for each plane
        all_faces: List of face arrays for each plane
    """
    try:
        if not all_vertices or not all_faces:
            print(f"Warning: No planes to write to {filepath}")
            return

        # Combine all vertices and faces
        combined_vertices = []
        combined_faces = []
        vertex_offset = 0

        for vertices, faces in zip(all_vertices, all_faces):
            combined_vertices.append(vertices)
            # Offset face indices to account for previous vertices
            combined_faces.append(faces + vertex_offset)
            vertex_offset += len(vertices)

        # Concatenate
        final_vertices = np.vstack(combined_vertices)
        final_faces = np.vstack(combined_faces)

        # Write combined mesh
        write_cross_section_stl(filepath, final_vertices, final_faces)

        print(f"Wrote combined planes STL: {Path(filepath).name} ({len(all_vertices)} planes)")

    except Exception as e:
        raise IOError(f"Failed to write combined STL file {filepath}: {str(e)}")


def save_results_pickle(filepath: str, results: Any) -> None:
    """
    Save results to pickle file (replaces .mat file saving)

    Args:
        filepath: Output pickle file path
        results: Results object to save
    """
    import pickle

    try:
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, 'wb') as f:
            pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"Saved results to: {Path(filepath).name}")

    except Exception as e:
        raise IOError(f"Failed to save pickle file {filepath}: {str(e)}")


def load_results_pickle(filepath: str) -> Any:
    """
    Load results from pickle file

    Args:
        filepath: Input pickle file path

    Returns:
        Loaded results object
    """
    import pickle

    try:
        with open(filepath, 'rb') as f:
            results = pickle.load(f)

        print(f"Loaded results from: {Path(filepath).name}")
        return results

    except Exception as e:
        raise IOError(f"Failed to load pickle file {filepath}: {str(e)}")
