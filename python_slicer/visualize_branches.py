#!/usr/bin/env python3
"""Quick visualization of branch detection results."""

import sys
from pathlib import Path
import pyvista as pv

pv.OFF_SCREEN = True

output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./branch_output")
surface_stl = sys.argv[2] if len(sys.argv) > 2 else None

colors = {
    "LeftNose": "blue",
    "RightNose": "red",
    "Mouth": "green",
    "DescendingAirway": "gold",
}

# Load centerlines
centerlines = {}
for name, color in colors.items():
    cl_files = list((output_dir / name).glob("*_centerline.vtp"))
    if cl_files:
        centerlines[name] = pv.read(str(cl_files[0]))

# Load surface meshes
meshes = {}
for name, color in colors.items():
    stl_files = list((output_dir / name).glob("*_mesh.stl"))
    if stl_files:
        meshes[name] = pv.read(str(stl_files[0]))

# View 1: Centerlines on transparent original surface
p1 = pv.Plotter(off_screen=True, window_size=[1600, 1000])
if surface_stl:
    p1.add_mesh(pv.read(surface_stl), color="lightgray", opacity=0.12)
for name, cl in centerlines.items():
    p1.add_mesh(cl.tube(radius=0.8), color=colors[name], label=name)
p1.add_legend(size=(0.2, 0.2))
p1.camera_position = "yz"
p1.camera.zoom(1.2)
p1.screenshot(str(output_dir / "viz_centerlines_side.png"))
p1.camera_position = "xz"
p1.camera.zoom(1.2)
p1.screenshot(str(output_dir / "viz_centerlines_front.png"))
p1.close()

# View 2: Split surface meshes
p2 = pv.Plotter(off_screen=True, window_size=[1600, 1000])
for name, mesh in meshes.items():
    p2.add_mesh(mesh, color=colors[name], opacity=0.6, label=name)
p2.add_legend(size=(0.2, 0.2))
p2.camera_position = "yz"
p2.camera.zoom(1.2)
p2.screenshot(str(output_dir / "viz_surface_side.png"))
p2.camera_position = "xz"
p2.camera.zoom(1.2)
p2.screenshot(str(output_dir / "viz_surface_front.png"))
p2.close()

print(f"Saved visualizations to {output_dir}/viz_*.png")
