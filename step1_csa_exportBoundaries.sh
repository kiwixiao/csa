#!/bin/bash
source /storage/Qiwei/shellScripts/check_status.sh


echo "This script runs Starccm+ in batch mode"
echo "Please make sure your final simulation saved as a name with 'final' in it, so the script\n will locate it precisely without issue."
echo "This script is adjusted for Inpsire project, maybe not working directly for other projects."
read -e -p "what is the subject ID: e.g. Inspire02 or Inspire03 etc: " sub
#read -e -p "First macro to run by default is exporting boundaries from simulation file\n: (input is optional, it has default value) " macro
#read -e -p "The second macro need to run is exploting star data plots from simulation file to allow\n analysis of simulation outputs." macro2

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


macro="${macro:=./functions/exportBoundaryForCSA.java}"
macro2="${macro2:=./functions/ExportStarPlotData.java}"

#		replace the subject name for: export boundary macro
sed -i "s/private\ String\ sub\ = .*;/private\ String\ sub\ =\ \"$sub\";/" $macro
# 		replace the category string using the input data: export boundary macro
sed -i "s/private\ String\ category\ = .*;/private\ String\ category\ =\ \"$category\";/" $macro

#		replace the subject name with correct subject name in: ExportStarPlotData macro
sed -i "s/private\ String\ subject\ = .*;/private\ String\ subject\ =\ \"$sub\";/" $macro2
#		replace the category string using input string: ExportStarPlotData.java macro
sed -i "s/private\ String\ category\ = .*;/private\ String\ category\ =\ \"$category\";/" $macro2

#		define the simulation path
simpath="/storage/Qiwei/cchmc_OSA/Projects_local/InspireProject/$sub$category/$sub*final*.sim"

echo This is the macro being picked: $macro
echo This is the macro being picked second: $macro2

/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro -power -podkey yDmwnSTsDqC1mE54Rl4TIw -licpath 1999@flex.cd-adapco.com -rsh ssh $simpath

./functions/convertSTLwithScaledPostfix.sh

/opt/Siemens/14.06.012-R8/STAR-CCM+14.06.012-R8/star/bin/starccm+ -batch $macro2 -power -podkey yDmwnSTsDqC1mE54Rl4TIw -licpath 1999@flex.cd-adapco.com -rsh ssh $simpath
