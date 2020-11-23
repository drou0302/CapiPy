#!/bin/bash
echo "Indicate the full location of conda.sh"
read $CONDA 
source $CONDA

echo Creating vCapPy environment...
conda env create -f ./Environment/vCapiPy.yml > ./env_creation.log 2>&1
echo Environment created.

echo Copying files...
cp Environment/*.sh ../
echo "Enter your MODELLER license (https://salilab.org/modeller/registration.html):"
read $LICENSE 
python ./Environment/file_modification.py $1 $LICENSE >> copyfiles1.log 2>&1

cd Environment
touch installation.log
cd ..
cp ./Environment/installation.log ./
echo Finished copying files.

read -rsp $'Press enter to close the finish.\n'

