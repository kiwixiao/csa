#!/bin/bash

echo "-- this code use mirtk covert vtp to vtk for later use"

# i assume this step can be avoided if we change the input file format of the matlab script.

echo "check and create vtk folder"

if [ -d "vtk" ]; then
	echo "vtk folder already exist"
	echo "deleting now"
	rm -rf "vtk"
	echo "deleted"
	echo "creating new one"
	mkdir "vtk"
	echo 'done'
else
	echo "vtk folder not exist"
	echo "creating now"
	mkdir "vtk"
	echo 'done'
fi

for file in ./vtp/*.vtp;
do
	echo "loop through all the vtp convert to vtk for matlab"
	echo "converting file "$file" "
	mirtk transform-points $file "${file%.*}.vtk" -ascii;
	echo "moving vtk file to vtk folder"
	mv ./vtp/*.vtk ./vtk/
done
	
