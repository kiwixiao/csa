#!/bin/bash
echo "This is calling matlab code to calculate the diameter of CSA planes"
read -e -p "Which side of the geometry is calculating: if left input 1, \nif right side input 2" side
#read -e -p 'what is the subject file folder name: such as:(LeftNoseDecendingSlicedSTLs or RightNoseSlicedSTLs): ' filePath
#read -e -p 'what is the number of timePoints interpolated from using mirtk: ' noOfTime

matlab -nodisplay -r " try cd('../../slicerMatlabFunctions'); PerimeterandhydDiamCalc("$side"); catch; end; quit"
