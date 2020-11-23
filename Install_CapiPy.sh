#!/bin/bash
echo "Indicate the full location of conda.sh"
read $CONDA 
source $CONDA
echo Creating vCapPy environment...
conda env create -f ./Environment/vCapiPy.yml > ./ENV-CREATION.log 
echo Environment created.
echo ____________________________________________________________________________

cp Environment/*.sh ../
echo "Enter your MODELLER license (https://salilab.org/modeller/registration.html):"
read $LICENSE 
python ./Environment/file_modification.py $1 $LICENSE > CON
echo File modification finished...
echo ____________________________________________________________________________

read -rsp $'Press enter to close the finish.\n'

