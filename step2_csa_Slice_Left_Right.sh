#!/bin/bash

echo this script calling scrpts in folder LeftNoseDecending/FFD/

# the round bracket will create a subshell, so the pwd will be the same as cd return, 
# after it is done, it also does not affect the pwd in this main shell.
(cd ./LeftNoseDecending/FFD && ./step2To5.sh)
if [ $? eq 0 ]; then
	echo "Left side slicing success."
else
	echo "Left side slicing failed."
	echo "Existing the code"
	exit 1
fi

(cd ./RightNose/FFD && ./step2To5.sh)
if [ $? eq 0 ]; then
	echo Right slicing code done successfully.
else
	echo "Right slicing code failed."
	exit 1
fi
