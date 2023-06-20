#!/bin/bash
#source /storage/Qiwei/shellScripts/check_status.sh

check_status(){
	# this function takes the first positional argument as the error message,
	# it will move nex line if there is no error, else it will exist the code
	if [ $? -ne 0 ]; then
		echo $1
		exit 1
	else
		echo -e "\[31m[SUCCESS]\e[0m"
	fi	
}

echo "This script runs Starccm+ in batch mode"
echo "Please make sure use the final simulation if you want to export the star plot table as well for non-geometric analysis purpose"

echo "This script is adjusted for Inpsire project, maybe not working directly for other projects."
read -e -p "what is the subject ID: e.g. Inspire02 or Inspire03 etc: " sub

if echo $sub | grep -q "Inspire"; then
	echo "This is for Inspire project."
	read -e -p "what is the image number of this subject: e.g.: CT003 or CT006 etc: " category
elif echo $sub | grep -q "OSAMRI"; then
	echo "This is for OSAMRI project"
	echo "Will set the category variable to empty for compatibility"
	category=""
else
	echo wrong input for subject name.
	echo existing code with code 11
	exit 11
fi

read -e -p "Where is the simulation you have to use, please give the full path: " simpath

# also please make sure this shell script is in the same folder with another folder: functions
check_dir(){
	# check is a dir exists or not
	if [ -d "$1" ]; then
		echo "Directory $1 exists."
		echo "this is essential if you do not want to manually worry about\n where are the java macros"
	else
		echo "Directory $1 does not in the same location as the current shell script"
		echo "Stop the code now"
		exit 1
	fi
}
check_dir "./functions"
check_status "functions folder not in the same location as the current shell script"

read -e -p "Do you want to export 3D boundaries and the star plot table at the same time? (yes/no): " ans

macro="${macro:=./functions/exportBoundaryForCSA.java}"
macro2="${macro2:=./functions/ExportStarPlotData.java}"

if [ $ans == "No|no|N|n" ]; then
	#	replace the subject name for: export boundary macro
	sed -i "s/private\ String\ sub\ = .*;/private\ String\ sub\ =\ \"$sub\";/" $macro
	check_status "replacing subject ID in export boundary macro failed."
	# 	replace the category string using the input data: export boundary macro
	sed -i "s/private\ String\ category\ = .*;/private\ String\ category\ =\ \"$category\";/" $macro

else
	#	replace the subject name for: export boundary macro
	sed -i "s/private\ String\ sub\ = .*;/private\ String\ sub\ =\ \"$sub\";/" $macro
	check_status "replacing subject ID in export boundary macro failed."
	# 	replace the category string using the input data: export boundary macro
	sed -i "s/private\ String\ category\ = .*;/private\ String\ category\ =\ \"$category\";/" $macro

#	replace the subject name with correct subject name in: ExportStarPlotData macro
	sed -i "s/private\ String\ subject\ = .*;/private\ String\ subject\ =\ \"$sub\";/" $macro2
#	replace the category string using input string: ExportStarPlotData.java macro
	sed -i "s/private\ String\ category\ = .*;/private\ String\ category\ =\ \"$category\";/" $macro2
fi

#		define the simulation path
#simpath="/storage/Qiwei/cchmc_OSA/Projects_local/InspireProject/$sub$category/$sub*final*.sim"

echo This is marco for export boundaries from CFD simulation: $macro
echo This is the macro for export star plot table for non-geometric analysis: $macro2

locate_file(){
	if [ ! -f "$1" ]; then
		echo "File $1 not found"
		exit 1
	else
		echo "Success: Found file $1"
		echo get the file content
		lic_key=$(cat $1)
	fi
}

locate_file "/home/xiaz9n/star_license_key.txt"

/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro -power -podkey $lic_key -licpath 1999@flex.cd-adapco.com -rsh ssh $simpath
check_status "Exporting 3D geometries failed, it can be either wrong macro setup, or it can bestar license issues, look at the output carefully"

./functions/convertSTLwithScaledPostfix.sh
check_status "Scale stl back to meters failed"

/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro2 -power -podkey $lic_key -licpath 1999@flex.cd-adapco.com -rsh ssh $simpath
