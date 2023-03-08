#!/bin/bash

echo this script calling scrpts in folder LeftNoseDecending/FFD/

(cd ./LeftNoseDecending/FFD && ./step2To4_unified.sh)

(cd ./RightNose/FFD && ./step2To4_unified.sh)
