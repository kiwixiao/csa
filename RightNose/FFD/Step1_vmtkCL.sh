#!/bin/bash
read -e -p "what is the input stl name: " instl
#read -e -p "what is the reducing factor: " f

cl="${instl%.*}_CL.vtp"

vmtkcenterlines -resampling 1 -resamplingstep 2 -ifile $instl -ofile $cl
#vmtkcenterlines -seedselector openprofiles -endpoints 1 -ifile $instl -resampling 1 -resamplingstep 2 -ofile $cl

vmtkcenterlinegeometry -ifile $cl -smoothing 1 -iterations 30 -factor 0.2 -outputsmoothed 1 -ofile "${cl%.*}_smooth.dat"
vmtkcenterlinegeometry -ifile $cl -smoothing 1 -iterations 30 -factor 0.2 -outputsmoothed 1 -ofile "${cl%.*}_smooth.vtp"

mv $cl "${cl%.*}"_original.vtp

# now we are running mesh_slicer to get the plane sliced.
# first of all get the centerline coordinates and normals
cut -d " " -f 1,2,3,10,11,12 "${cl%.*}_smooth.dat" > "${cl%.*}_corNor.dat"

# line below designed to work with both mac and unix
lines=$(wc -l "${cl%.*}_corNor.dat" | grep -o -E '[0-9]+' | head -1 | sed -e 's/^0\+//')

echo " number of lines of files $lines"
sed "1s/.*/"$lines"/" "${cl%.*}_corNor.dat" > "${cl%.*}_Planes.dat"

#sed -n -e "1~"$f"p" "${cl%.*}_Planes.dat" > "${cl%.*}_Planes_Reduced.dat"
