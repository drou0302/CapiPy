#!/bin/bash
echo Creating vCapPy environment...
conda env create -f vCapiPy.yml > ./Environment/env_creation.log 2>&1
echo Environment created.
echo Copying files...
cp Environment/shfiles/*.sh ../
set /p LICENSE = "Enter your MODELLER license ((https://salilab.org/modeller/registration.html):"
python ./Environment/file_modification.py  %1 %LICENSE% >> ./Environment/copyfiles2.log 2>&1
cd Environment
touch installation.log
copy *.log installation.log
cd ..
cp ./Environment/installation.log ./
echo Finished copying files.
pause