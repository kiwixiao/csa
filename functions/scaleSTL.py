# import lib

import trimesh
import argparse
import os

# create argument parser
parser = argparse.ArgumentParser()
# add a positional arguments
parser.add_argument("input_file",help="The input file to be processed")
#parser the argument
args = parser.parse_args()
# now you can use input_file argument in your script
print("Processing file",args.input_file)
# get the input mesh path and file
path = args.input_file
# extract the file name from the input
names = path.split('/') # this gives all the substring separated by /
file_name = names[-1] # take the last index which is the filename
print("This is the file name from the input_file", file_name)

# load the stl
mesh = trimesh.load_mesh(path)
# creating a scaling factor
scale_factor = [1000,1000,1000]
# creating a scaleing transform using scale_factor
scaler = trimesh.transformations.scale_matrix(scale_factor)
# apply the transform to the mesh
mesh_scaled = mesh.apply_scale(1000)
# export the scaled mesh.
# append scaled to the original mesh name
out_file_name = file_name.split('.')[0]+'_Scaled.stl'
mesh_scaled.export(out_file_name)

