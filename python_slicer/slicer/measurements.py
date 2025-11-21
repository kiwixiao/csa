"""
Geometric measurements for cross-sections
Replaces: DiameterCalc.m, PerimeterandhydDiamCalc.m
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Optional
from dataclasses import dataclass
from shapely.geometry import Polygon
from scipy.optimize import least_squares


@dataclass
class CrossSectionMetrics:
    """
    Container for all measurements of a single cross-section
    """
    arc_length: float  # Position along centerline (mm)
    area: float  # Cross-sectional area (mm^2)
    perimeter: float  # Perimeter length (mm)
    centroid: np.ndarray  # 3D centroid position
    hydraulic_diameter: float  # 4 * Area / Perimeter
    equivalent_diameter: float  # Diameter of circle with same area
    major_axis: float  # Ellipse major axis
    minor_axis: float  # Ellipse minor axis
    diameter_ratio: float  # Major / Minor axis (elongation)
    plane_index: int  # Index in original centerline
    is_valid: bool  # Quality flag


class DiameterCalculator:
    """
    Calculate various diameter metrics for cross-sections

    Replaces: DiameterCalc.m and PerimeterandhydDiamCalc.m
    """

    @staticmethod
    def compute_hydraulic_diameter(area: float, perimeter: float) -> float:
        """
        Compute hydraulic diameter: Dh = 4 * Area / Perimeter

        This is commonly used in fluid dynamics to characterize non-circular ducts

        Args:
            area: Cross-sectional area (mm^2)
            perimeter: Perimeter (mm)

        Returns:
            Hydraulic diameter (mm)
        """
        if perimeter <= 0:
            return 0.0
        return 4.0 * area / perimeter

    @staticmethod
    def compute_equivalent_diameter(area: float) -> float:
        """
        Compute diameter of circle with equivalent area: D = sqrt(4*A/pi)

        Args:
            area: Cross-sectional area (mm^2)

        Returns:
            Equivalent circular diameter (mm)
        """
        if area <= 0:
            return 0.0
        return np.sqrt(4.0 * area / np.pi)

    @staticmethod
    def fit_ellipse_to_points(points: np.ndarray) -> Tuple[float, float, np.ndarray, float]:
        """
        Fit an ellipse to 2D boundary points and extract parameters

        Args:
            points: Nx2 array of boundary points

        Returns:
            Tuple of (major_axis, minor_axis, center, rotation_angle)
        """
        try:
            # Fit ellipse using least squares
            # Ellipse equation: ((x-cx)/a)^2 + ((y-cy)/b)^2 = 1
            # After rotation by theta

            def ellipse_residuals(params, points):
                cx, cy, a, b, theta = params
                cos_t = np.cos(theta)
                sin_t = np.sin(theta)

                # Rotate points to ellipse frame
                xr = (points[:, 0] - cx) * cos_t + (points[:, 1] - cy) * sin_t
                yr = -(points[:, 0] - cx) * sin_t + (points[:, 1] - cy) * cos_t

                # Residuals from ellipse equation
                residuals = (xr / a) ** 2 + (yr / b) ** 2 - 1.0
                return residuals

            # Initial guess
            cx0 = np.mean(points[:, 0])
            cy0 = np.mean(points[:, 1])
            a0 = np.std(points[:, 0])
            b0 = np.std(points[:, 1])
            theta0 = 0.0

            initial = [cx0, cy0, a0, b0, theta0]

            # Fit
            result = least_squares(ellipse_residuals, initial, args=(points,))

            if result.success:
                cx, cy, a, b, theta = result.x

                # Ensure a >= b (a is major axis)
                if a < b:
                    a, b = b, a
                    theta += np.pi / 2

                major_axis = 2 * abs(a)
                minor_axis = 2 * abs(b)
                center = np.array([cx, cy])

                return major_axis, minor_axis, center, theta
            else:
                # Fallback: use PCA
                return DiameterCalculator._fit_ellipse_pca(points)

        except Exception:
            # Fallback to PCA method
            return DiameterCalculator._fit_ellipse_pca(points)

    @staticmethod
    def _fit_ellipse_pca(points: np.ndarray) -> Tuple[float, float, np.ndarray, float]:
        """
        Fit ellipse using PCA (Principal Component Analysis)

        Args:
            points: Nx2 array of boundary points

        Returns:
            Tuple of (major_axis, minor_axis, center, rotation_angle)
        """
        # Center points
        center = np.mean(points, axis=0)
        centered_points = points - center

        # Compute covariance matrix
        cov = np.cov(centered_points.T)

        # Eigenvalue decomposition
        eigenvalues, eigenvectors = np.linalg.eig(cov)

        # Sort by eigenvalues (descending)
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Major and minor axes (2 * std along principal components)
        major_axis = 2 * np.sqrt(eigenvalues[0])
        minor_axis = 2 * np.sqrt(eigenvalues[1])

        # Rotation angle
        theta = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])

        return major_axis, minor_axis, center, theta

    @staticmethod
    def compute_diameter_ratio(major_axis: float, minor_axis: float) -> float:
        """
        Compute diameter ratio (elongation): ratio = major / minor

        This quantifies how elongated the cross-section is
        - Ratio = 1.0: perfectly circular
        - Ratio > 1.0: elongated (elliptical)

        Args:
            major_axis: Major axis length
            minor_axis: Minor axis length

        Returns:
            Diameter ratio
        """
        if minor_axis <= 0:
            return 0.0
        return major_axis / minor_axis


class DiameterProfile:
    """
    Analyze diameter profile along entire centerline

    Replaces: Combined logic from DiameterCalc.m and data assembly in SlicerMasterCode.m
    """

    def __init__(self):
        """Initialize empty profile"""
        self.metrics_list: List[CrossSectionMetrics] = []

    def add_cross_section(self,
                         arc_length: float,
                         area: float,
                         perimeter: float,
                         centroid: np.ndarray,
                         boundary_2d: np.ndarray,
                         plane_index: int,
                         is_valid: bool = True) -> None:
        """
        Add a cross-section to the profile

        Args:
            arc_length: Position along centerline
            area: Cross-sectional area
            perimeter: Perimeter length
            centroid: 3D centroid
            boundary_2d: 2D boundary points for ellipse fitting
            plane_index: Index in centerline
            is_valid: Quality flag
        """
        # Compute basic diameters
        hydraulic_diam = DiameterCalculator.compute_hydraulic_diameter(area, perimeter)
        equivalent_diam = DiameterCalculator.compute_equivalent_diameter(area)

        # Fit ellipse and compute axes
        try:
            major_axis, minor_axis, _, _ = DiameterCalculator.fit_ellipse_to_points(boundary_2d)
            diameter_ratio = DiameterCalculator.compute_diameter_ratio(major_axis, minor_axis)
        except Exception:
            # Fallback: use equivalent diameter
            major_axis = equivalent_diam
            minor_axis = equivalent_diam
            diameter_ratio = 1.0

        # Create metrics object
        metrics = CrossSectionMetrics(
            arc_length=arc_length,
            area=area,
            perimeter=perimeter,
            centroid=centroid.copy(),
            hydraulic_diameter=hydraulic_diam,
            equivalent_diameter=equivalent_diam,
            major_axis=major_axis,
            minor_axis=minor_axis,
            diameter_ratio=diameter_ratio,
            plane_index=plane_index,
            is_valid=is_valid
        )

        self.metrics_list.append(metrics)

    def to_dataframe(self, flip_order: bool = False) -> pd.DataFrame:
        """
        Convert profile to pandas DataFrame

        Args:
            flip_order: If True, reverse order (for MATLAB compatibility)

        Returns:
            DataFrame with all measurements
        """
        if not self.metrics_list:
            return pd.DataFrame()

        # Build data dictionary
        data = {
            'plane_index': [m.plane_index for m in self.metrics_list],
            'arc_length_mm': [m.arc_length for m in self.metrics_list],
            'area_mm2': [m.area for m in self.metrics_list],
            'perimeter_mm': [m.perimeter for m in self.metrics_list],
            'hydraulic_diameter_mm': [m.hydraulic_diameter for m in self.metrics_list],
            'equivalent_diameter_mm': [m.equivalent_diameter for m in self.metrics_list],
            'major_axis_mm': [m.major_axis for m in self.metrics_list],
            'minor_axis_mm': [m.minor_axis for m in self.metrics_list],
            'diameter_ratio': [m.diameter_ratio for m in self.metrics_list],
            'centroid_x': [m.centroid[0] for m in self.metrics_list],
            'centroid_y': [m.centroid[1] for m in self.metrics_list],
            'centroid_z': [m.centroid[2] for m in self.metrics_list],
            'is_valid': [m.is_valid for m in self.metrics_list]
        }

        df = pd.DataFrame(data)

        if flip_order:
            df = df.iloc[::-1].reset_index(drop=True)

        return df

    def get_summary_statistics(self) -> dict:
        """
        Compute summary statistics across all cross-sections

        Returns:
            Dictionary with min, max, mean, std for key metrics
        """
        if not self.metrics_list:
            return {}

        areas = np.array([m.area for m in self.metrics_list if m.is_valid])
        hyd_diams = np.array([m.hydraulic_diameter for m in self.metrics_list if m.is_valid])
        ratios = np.array([m.diameter_ratio for m in self.metrics_list if m.is_valid])

        summary = {
            'area_min_mm2': np.min(areas) if len(areas) > 0 else 0,
            'area_max_mm2': np.max(areas) if len(areas) > 0 else 0,
            'area_mean_mm2': np.mean(areas) if len(areas) > 0 else 0,
            'area_std_mm2': np.std(areas) if len(areas) > 0 else 0,
            'hydraulic_diameter_min_mm': np.min(hyd_diams) if len(hyd_diams) > 0 else 0,
            'hydraulic_diameter_max_mm': np.max(hyd_diams) if len(hyd_diams) > 0 else 0,
            'hydraulic_diameter_mean_mm': np.mean(hyd_diams) if len(hyd_diams) > 0 else 0,
            'diameter_ratio_min': np.min(ratios) if len(ratios) > 0 else 0,
            'diameter_ratio_max': np.max(ratios) if len(ratios) > 0 else 0,
            'diameter_ratio_mean': np.mean(ratios) if len(ratios) > 0 else 0,
            'n_valid_planes': len([m for m in self.metrics_list if m.is_valid]),
            'n_total_planes': len(self.metrics_list)
        }

        return summary

    def find_minimum_area_location(self) -> Optional[CrossSectionMetrics]:
        """
        Find location of minimum cross-sectional area

        Useful for identifying constrictions

        Returns:
            CrossSectionMetrics at minimum area, or None if no valid sections
        """
        valid_metrics = [m for m in self.metrics_list if m.is_valid]

        if not valid_metrics:
            return None

        min_metric = min(valid_metrics, key=lambda m: m.area)
        return min_metric
