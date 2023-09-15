#!/bin/bash

echo 'This is calling matlab slicer code: SlicerMasterCode'

read -e -p 'what is the subject name: ' sub
read -e -p 'what the division path: such as: LeftNoseDecending or RightNose: ' path

matlab -nodisplay -r "try; cd('../../slicerMatlabFunctions'); SlicerMasterCode('$sub','$path'); catch ME; disp(ME.message); end; quit"
