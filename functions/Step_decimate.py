#!/usr/bin/python3
"""
Decimation
~~~~~~~~~~

Decimate a mesh

"""
# sphinx_gallery_thumbnail_number = 4
import subprocess
import sys
#import pyvista as pv
#from pyvista import examples
#import pyinputplus as pyip
import os
from stl import mesh
import numpy as np
import trimesh

# userin = pyip.inputStr(prompt='Enter the path and name of the stl File:')
userin = sys.argv[1]
print("this is the file for decimate")
#mesh = pv.read(userin)
mesh = trimesh.load_mesh(userin)
mesh.fix_normals()
mesh.decimate(target_percent=0.3)
mesh.export(userin[:-4]+"_decimated.stl")
#decimated.plot(cpos=cpos, **dargs)

