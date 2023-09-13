#!/bin/bash
echo 'this script is about to run through Step2 to Step 4, finger crosses this will work' > debug_masterShellCodeLogFile.txt

#read -e -p "subject name:(optional) " sub
#read -e -p "partition name: LeftNoseDecending or RightNose(optional): " path

sub=$(basename $(ls *Scaled.stl) | cut -d '_' -f 1)
path=$(basename $(ls *Scaled.stl) | cut -d '_' -f 2)

echo this is the subject name: $sub
echo this is the path for runMatlabSlicer input $path

./Step2_geoFromDof.sh && echo 'Step 2 done so far' >> debug_masterShellCodeLogFile.txt


./Step3_vtpTovtk.sh && echo 'Step 3 done so far' >> debug_masterShellCodeLogFile.txt


../../functions/runMatlabSlicer.sh "$sub" $path && echo 'Step 4 done now. this should be all if you see this.' >> debug_masterShellCodeLogFile.txt

./Step5_cal_Diameter.sh && echo Seeing this line means step 5 is done, FYI, the diameter calculation is also using default 100 ms as the interloplation step.
