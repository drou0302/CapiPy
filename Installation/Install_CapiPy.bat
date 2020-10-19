@echo off
echo Creating vCapiPy environment...
call conda env create -f Environment/vCapiPy.yml >> ./Environment/env_creation.log 2>&1
echo Environment created

echo Copying files...
copy ".\Environment\batFiles\CapiPy.bat" "..\" >> ./Environment/copyfiles.log 2>&1
set /p LICENSE="Enter your MODELLER license (https://salilab.org/modeller/registration.html):"
call python ./Environment/file_modification.py %1 %LICENSE% >> ./Environment/copyfiles2.log 2>&1
echo Finished copying files.

cd Environment
type NUL > installation.log
copy *.log  installation.log
copy installation.log ..\
echo Installation complete. Check installation.log for errors.
pause
