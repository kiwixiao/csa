"""
Tests for branch_detector module.

Unit tests for graph building, endpoint classification, and path tracing.
Integration tests with real STL require the vmtk conda environment and
user-provided data.
"""

import pytest
import numpy as np
import networkx as nx


# ---------------------------------------------------------------------------
# Test helpers - mock node_info structures
# ---------------------------------------------------------------------------

def make_node_info(gid, centroid, mean_radius, blanking=0):
    """Helper to create a mock node_info entry"""
    return {
        "points": np.array([centroid]),
        "radii": np.array([mean_radius]),
        "blanking": blanking,
        "centerline_ids": set(),
        "tract_ids": set(),
        "start": np.array(centroid),
        "end": np.array(centroid),
        "centroid": np.array(centroid),
        "mean_radius": mean_radius,
    }


# ---------------------------------------------------------------------------
# Test classify_endpoints
# ---------------------------------------------------------------------------

class TestClassifyEndpoints:
    """Test anatomical labeling of terminal branch segments"""

    def _import_classify(self):
        """Import classify_endpoints - skip if vmtk not available"""
        try:
            from branch_detector import classify_endpoints
            return classify_endpoints
        except ImportError:
            pytest.skip("vmtk not installed, skipping branch_detector tests")

    def test_three_terminals_two_nostrils_one_outlet(self):
        """Standard case: 2 nostrils + 1 outlet"""
        classify_endpoints = self._import_classify()

        G = nx.Graph()
        node_info = {}

        # Outlet: large radius, low Z (caudal)
        node_info[0] = make_node_info(0, [0, 0, -50], mean_radius=10.0)
        # Left nostril: small radius, high Z, negative X
        node_info[1] = make_node_info(1, [-15, 10, 50], mean_radius=3.0)
        # Right nostril: small radius, high Z, positive X
        node_info[2] = make_node_info(2, [15, 10, 50], mean_radius=3.0)

        for gid in [0, 1, 2]:
            G.add_node(gid, **node_info[gid])

        terminals = [0, 1, 2]
        labels = classify_endpoints(G, node_info, terminals)

        assert labels["DescendingAirway"] == 0
        assert labels["LeftNose"] == 1
        assert labels["RightNose"] == 2

    def test_four_terminals_two_nostrils_mouth_outlet(self):
        """Case with mouth: 2 nostrils + mouth + outlet"""
        classify_endpoints = self._import_classify()

        G = nx.Graph()
        node_info = {}

        # Outlet: large radius, lowest Z
        node_info[0] = make_node_info(0, [0, 0, -50], mean_radius=10.0)
        # Left nostril: high Z
        node_info[1] = make_node_info(1, [-15, 10, 50], mean_radius=3.0)
        # Right nostril: high Z
        node_info[2] = make_node_info(2, [15, 10, 50], mean_radius=3.0)
        # Mouth: medium Z, anterior
        node_info[3] = make_node_info(3, [0, 20, 10], mean_radius=4.0)

        for gid in [0, 1, 2, 3]:
            G.add_node(gid, **node_info[gid])

        terminals = [0, 1, 2, 3]
        labels = classify_endpoints(G, node_info, terminals)

        assert labels["DescendingAirway"] == 0
        assert labels["LeftNose"] == 1
        assert labels["RightNose"] == 2
        assert labels["Mouth"] == 3

    def test_outlet_identified_by_largest_radius(self):
        """Outlet should be the terminal with largest mean radius"""
        classify_endpoints = self._import_classify()

        G = nx.Graph()
        node_info = {}

        # Even if outlet is at high Z, it should be identified by radius
        node_info[0] = make_node_info(0, [0, 0, 100], mean_radius=15.0)  # largest radius
        node_info[1] = make_node_info(1, [-15, 10, 50], mean_radius=2.0)
        node_info[2] = make_node_info(2, [15, 10, 50], mean_radius=2.0)

        for gid in [0, 1, 2]:
            G.add_node(gid, **node_info[gid])

        terminals = [0, 1, 2]
        labels = classify_endpoints(G, node_info, terminals)

        assert labels["DescendingAirway"] == 0

    def test_left_right_distinguished_by_x(self):
        """Left nostril should have smaller X coordinate"""
        classify_endpoints = self._import_classify()

        G = nx.Graph()
        node_info = {}

        node_info[0] = make_node_info(0, [0, 0, -50], mean_radius=10.0)
        node_info[1] = make_node_info(1, [-20, 0, 50], mean_radius=3.0)
        node_info[2] = make_node_info(2, [20, 0, 50], mean_radius=3.0)

        for gid in [0, 1, 2]:
            G.add_node(gid, **node_info[gid])

        terminals = [0, 1, 2]
        labels = classify_endpoints(G, node_info, terminals)

        assert labels["LeftNose"] == 1   # x=-20 (smaller)
        assert labels["RightNose"] == 2  # x=20 (larger)


# ---------------------------------------------------------------------------
# Test trace_partition_paths
# ---------------------------------------------------------------------------

class TestTracePartitionPaths:
    """Test path tracing from terminals to outlet"""

    def _import_trace(self):
        try:
            from branch_detector import trace_partition_paths
            return trace_partition_paths
        except ImportError:
            pytest.skip("vmtk not installed, skipping branch_detector tests")

    def test_simple_y_topology(self):
        """Y-shaped airway: two nostrils merge into one descending"""
        trace_partition_paths = self._import_trace()

        # Graph: 1 -- 3 -- 0
        #         2 /
        G = nx.Graph()
        node_info = {}
        for gid in [0, 1, 2, 3]:
            node_info[gid] = make_node_info(gid, [0, 0, 0], 5.0)

        G.add_node(0)  # outlet
        G.add_node(1)  # left nose
        G.add_node(2)  # right nose
        G.add_node(3)  # shared descending segment
        G.add_edge(1, 3)
        G.add_edge(2, 3)
        G.add_edge(3, 0)

        labels = {
            "DescendingAirway": 0,
            "LeftNose": 1,
            "RightNose": 2,
        }

        partitions = trace_partition_paths(G, node_info, labels)

        # LeftNoseDecending should go 1 -> 3 -> 0
        assert partitions["LeftNoseDecending"] == [1, 3, 0]
        # RightNoseDecending should go 2 -> 3 -> 0
        assert partitions["RightNoseDecending"] == [2, 3, 0]
        # Short partitions
        assert partitions["LeftNose"] == [1]
        assert partitions["RightNose"] == [2]
        assert partitions["DescendingAirway"] == [0]

    def test_with_mouth_branch(self):
        """Y+T topology: two nostrils + mouth merge into descending"""
        trace_partition_paths = self._import_trace()

        # Graph: 1 -- 4 -- 5 -- 0
        #         2 /    3 /
        G = nx.Graph()
        node_info = {}
        for gid in [0, 1, 2, 3, 4, 5]:
            node_info[gid] = make_node_info(gid, [0, 0, 0], 5.0)

        for gid in [0, 1, 2, 3, 4, 5]:
            G.add_node(gid)
        G.add_edge(1, 4)  # left nose to first junction
        G.add_edge(2, 4)  # right nose to first junction
        G.add_edge(4, 5)  # first junction to second junction
        G.add_edge(3, 5)  # mouth to second junction
        G.add_edge(5, 0)  # second junction to outlet

        labels = {
            "DescendingAirway": 0,
            "LeftNose": 1,
            "RightNose": 2,
            "Mouth": 3,
        }

        partitions = trace_partition_paths(G, node_info, labels)

        assert partitions["LeftNoseDecending"] == [1, 4, 5, 0]
        assert partitions["RightNoseDecending"] == [2, 4, 5, 0]
        assert partitions["MouthDecending"] == [3, 5, 0]


# ---------------------------------------------------------------------------
# Test build_branch_graph (requires vmtk structures, so lighter tests)
# ---------------------------------------------------------------------------

class TestCenterlineToVtp:
    """Test centerline point array to VTP conversion"""

    def _import_func(self):
        try:
            from branch_detector import centerline_to_vtp
            return centerline_to_vtp
        except ImportError:
            pytest.skip("vmtk not installed")

    def test_basic_polyline(self):
        centerline_to_vtp = self._import_func()

        points = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
        ], dtype=float)

        vtp = centerline_to_vtp(points)
        assert vtp.GetNumberOfPoints() == 3
        assert vtp.GetNumberOfLines() == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
