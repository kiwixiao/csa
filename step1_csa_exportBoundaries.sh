#!/bin/bash
source /storage/Qiwei/shellScripts/check_status.sh


echo "This script runs Starccm+ in batch mode"
#echo "Please make sure your final simulation saved as a name with 'final' in it, so the script\n will locate it precisely without issue."
echo "This script is adjusted for Inpsire project, maybe not working directly for other projects."
read -e -p "what is the subject ID: e.g. Inspire02 or OSAMRI016 etc: " sub
#read -e -p "First macro to run by default is exporting boundaries from simulation file\n: (input is optional, it has default value) " macro
#read -e -p "The second macro need to run is exploting star data plots from simulation file to allow\n analysis of simulation outputs." macro2
read -e -p "Simulation for stl export, which should be time zero simulation: " simpath
echo -e "Please enter the sim name for flow data export"
echo -e "Press enter to escape this step."
read -e -p "Simulation for flow table results, should be the final simulation: " simfinal
read -e -p "Is this a CPAP related simulation?[yes|no] " res
if [[ "$res" == "yes" ]]; then
	cpap="_CPAP"
else
	cpap=""
fi

if echo $sub | grep -q "Inspire"; then
	echo "This is for Inspire project."
	read -e -p "what is the image number of this subject: e.g.: CT003 or CT006 etc: " category
elif echo $sub | grep -q "OSAMRI"; then
	echo "This is for OSAMRI project"
	echo "Will set the category variable to empty for compatibility"
	category=""
else
	echo not correct input for subject name input.
	echo existing code with 1
	exit 1
fi
#folder for saving stl boundary files
stl_out_path="$(pwd)"
echo "This is the current working dir: "
echo "$stl_out_path"

# create the folder for flow results
echo "check if the folder for saving flow results exist"
folder="$stl_out_path/files_${sub}${cpap}" 
echo "This is the folder for flow results:"
echo "$folder"

if [ -d "$folder" ]; then
	echo "Folder $folder already exists, Not creating it."
else
	mkdir "$folder"
	echo "Folder $Folder created."
fi

macro="$stl_out_path/functions/exportBoundaryForCSA.java"
macro2="$stl_out_path/functions/ExportStarPlotData.java"

#		replace the subject name for: export boundary macro
sed -i "s#private String sub[ ]*=.*\;#private String sub = \"$sub\"\;#" "$macro"

sed -i "s#private String outPath[ ]*=.*\;#private String outPath = \"$stl_out_path\"\;#" "$macro"

sed -i "s#private String category[ ]*=.*\;#private String category = \"$category\"\;#" "$macro"


#		replace the subject name with correct subject name in: ExportStarPlotData macro
sed -i "s#private String subject = *.*\;#private String subject = \"$sub\"\;#" "$macro2"
#		replace the category string using input string: ExportStarPlotData.java macro
sed -i "s#private String category = *.*\;#private String category = \"$category\"\;#" "$macro2"
sed -i "s#private String outPath = *.*\;#private String outPath = \"$folder\/$sub$cpap\"\;#" "$macro2"

#		define the simulation path
#simpath="/storage/Qiwei/cchmc_OSA/Projects_local/InspireProject/$sub$category/$sub*final*.sim"

echo This is the macro being picked: $macro
echo This is the macro being picked second: $macro2

echo This is the real path of the simulation picked
abs_simpath=$(realpath "$simpath")
echo "$abs_simpath"
abs_simfinal=$(realpath "$simfinal")
echo "this is the simulation picked for flow table export:"
echo "$abs_simfinal"

/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro -power -podkey GRLLLPgPZLUaFnZBOtU8pw -licpath 1999@flex.cd-adapco.com -rsh ssh $abs_simpath

./functions/convertSTLwithScaledPostfix.sh

if [ -z "$simfinal" ]; then
	echo "No simulation input for flow data export, skipping this step!"
else
	echo "Exporting flow data to folder:"
	/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro2 -power -podkey GRLLLPgPZLUaFnZBOtU8pw -licpath 1999@flex.cd-adapco.com -rsh ssh $abs_simpath
fi

echo "Script finished"
