#!/bin/bash
echo "This is calling matlab code to calculate the diameter of CSA planes"

matlab -nodisplay -r "try; cd('./slicerMatlabFunctions'); PerimeterandhydDiamCalc(1); catch ME; disp(ME.message); end; exit"

if [ $? -eq 0 ]; then
	echo Left side perimeter calculation failed
	exit 1
fi

matlab -nodisplay -r "try cd('./slicerMatlabFunctions'); PerimeterandhydDiamCalc(2); catch ME; disp(ME.message); end; exit"

if [ $? -eq 0 ]; then
	echo Right side perimeter calculation failed
	exit 1
fi
