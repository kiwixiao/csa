#!/bin/bash
read -e -p "what is the input stl name: " instl

cl="${instl%.*}_CL.vtp"

vmtkcenterlines -ifile $instl -ofile $cl

vmtkcenterlinegeometry -ifile $cl -smoothing 1 -iterations 30 -factor 0.2 -outputsmoothed 1 -ofile "${cl%.*}_smooth.vtp"

rm $cl
