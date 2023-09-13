import sys
import pyvista as pv

# Get the input file name from the command line arguments
input_file = sys.argv[1]

# Load the STL file
#mesh = o3d.io.read_triangle_mesh(input_file)
mesh = pv.read(input_file)


# Decimate the mesh
decimated_mesh = mesh.decimate(target_reduction=0.95)  # Adjust the target_reduction value as needed

# Save the decimated mesh
output_file = input_file[:-4] + "_decimated.stl"

decimated_mesh.save(output_file)

