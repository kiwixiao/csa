"""
Python Slicer Package
Replaces MATLAB slicer pipeline for airway cross-sectional analysis
"""

__version__ = "1.0.0"

from .io_utils import (
    read_stl,
    read_vtk_centerline,
    write_cross_section_stl,
    write_measurements_csv
)

from .geometry import (
    compute_plane_normal,
    rotate_to_xy_plane,
    rotate_back_to_3d,
    compute_arc_length,
    validate_planes_no_intersection,
    check_adjacent_planes_intersect,
    check_plane_vertices_same_side
)

from .validation import (
    check_plane_enclosed_by_mesh,
    validate_centroid_proximity
)

from .plane_slicer import AirwaySlicer
from .measurements import CrossSectionMetrics, DiameterProfile

__all__ = [
    'read_stl',
    'read_vtk_centerline',
    'write_cross_section_stl',
    'write_measurements_csv',
    'compute_plane_normal',
    'rotate_to_xy_plane',
    'rotate_back_to_3d',
    'compute_arc_length',
    'validate_planes_no_intersection',
    'check_plane_enclosed_by_mesh',
    'validate_centroid_proximity',
    'check_adjacent_planes_intersect',
    'check_plane_vertices_same_side',
    'AirwaySlicer',
    'CrossSectionMetrics',
    'DiameterProfile'
]
