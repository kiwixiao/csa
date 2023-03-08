#!/bin/bash
read -e -p 'what is the input stl: ' stl

python3 ../../functions/Step_decimate.py "$stl"

# the scale function using meshlab server is not working anymore due to deperication by meshlab as well as x window forwarding between mac and ubuntu.

#meshlabserver -i *decimated.stl -o "${stl%.*}_decimated_meter.stl" -s ../../functions/meshlab_Scale1000.mlx
