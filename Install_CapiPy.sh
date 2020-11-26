#!/bin/bash
echo "Indicate the full location of conda.sh"
read $CONDA 
source $CONDA
echo Creating vCapPy environment...
conda env create -f ./ncfiles/vCapiPy.yml > ./ENV-CREATION.log 
echo Environment created.
echo ____________________________________________________________________________

mv ./ncfiles/CapiPy.sh ./
echo "Enter your MODELLER license (https://salilab.org/modeller/registration.html):"
read $LICENSE 
python ./ncfiles/file_modification.py $1 $LICENSE > CON
echo File modification finished...
echo ____________________________________________________________________________

read -rsp $'Press enter to close the finish.\n'

