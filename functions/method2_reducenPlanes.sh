#!/bin/bash

echo 'if this not working correctly, very likely due to sed on mac is not gun-sed.'
read -e -p "what is the reducing factor of nPlanes: " f

# line below is reduce the file by evey $f rows:
sed -n -e "1~"$f"p" *Planes.dat > F_"$f"_nPlanesInfo.dat

# now count how many lines of this file
#lines=$(wc -l F_"$f"_nPlanesInfo.dat | cut -d " " -f 6)
lines=$(wc -l < F_"$f"_nPlanesInfo.dat)

echo 'this is the number of lines:' $lines

# now change the first line if the file:
sed -i "1s/.*/$(( $lines-1 ))/g" F_"$f"_nPlanesInfo.dat
