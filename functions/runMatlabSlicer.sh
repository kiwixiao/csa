#!/bin/bash

echo 'This is calling matlab slicer code: SlicerMasterCode'

#read -e -p 'what is the subject name: ' sub
#read -e -p 'what the division path: such as: LeftNoseDecending or RightNose: ' path
echo "this will read the position argument from the master shell script."

#echo this is the working directory:
#filepath=$(pwd)
#echo $filepath

matlab -nodisplay -r " try cd('../../slicerMatlabFunctions'); SlicerMasterCode('$1','$2'); catch; end; quit"
