#!/bin/bash

echo "this script converts all the stl exprot from star to right units with Scaled post fix"

for f in *Scaled.stl
do
if [[ -f "${f}" ]]; then
	echo Scaled stl exist, will exist this code do nothing
	exit 1
fi
done

fileList=($(ls -d *.stl))
echo ${fileList[@]}

echo looping stl for python scale code
for stl in ${fileList[@]}
do
	echo $stl
	python ./functions/scaleSTL.py $stl
done


