# use dof file export stl use for Virtual DISE
#!/bin/bash
#read -e -p "Interpolation step: " inteStep
#read -e -p "Image at frame 0: " Image0
#read -e -p "Centerline frame 0: " centerline0
#read -e -p "STL frame 0: " stl0
#read -e -p "Do you have the input.txt file ready? " ans

#source ../../LeftNoseDecending/FFD/Step2_geoFromDof.sh 
#inteStep=$it

inteStep=100 # this is default value, consistant with left side.
Image0=$(ls *.nii.gz)
centerline0=$(ls *CL_smooth.vtp)
stl0=$(ls *Scaled.stl)

echo -e "\n\n Now the code is running"

# read the input.txt file to find out time info
	beginPoint=$(sed -n 2p input.txt | sed 's/.*,//')
	timePoints=$(sed -n 3p input.txt | sed 's/.*,//')
	breTime=$(sed -n 4p input.txt | sed 's/.*,//')
	dt=$(($breTime/$(($timePoints-1))));
	timeEnd=$(($dt*$(($timePoints-1))));
# firstly interpolate the centerline based
python3 ../../functions/centerlineInterp.py \
  --target "$Image0" \
  --dofs ffds.csv \
  --mesh "$centerline0" \
  --start 0 \
  --stop "$timeEnd" \
  --step "$inteStep" \
  --output-mesh "out_{t:06.0f}.vtp" \
  --output-table "motiontable_CTL.csv"

python3 ../../functions/Interp.py \
  --target "$Image" \
  --dofs ffds.csv \
  --mesh "$stl0" \
  --start 0 \
  --stop "$timeEnd" \
  --step "$inteStep" \
  --output-mesh "out_{t:06.0f}.stl" \
  --output-table "motiontable.csv"

# now make a directory to save the centerline
if [ -d "vtp" ]; then
	rm -rf "vtp"
	mkdir "vtp"
	echo "vtp folder exist, but re-created again"
else
	mkdir "vtp"
	echo "vtp folder do not exist before, but created a new one"
fi

if [ -d "stl" ]; then
	rm -rf "stl"
	mkdir "stl"
	echo "stl folder exist, but re-created again"
else
	mkdir "stl"
	echo "stl folder do not exist before, but created a new one"
fi

mv out*.vtp ./vtp
mv out*.stl ./stl

