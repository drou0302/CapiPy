@echo off
ECHO Step 1/3: Creating vCapiPy environment...
CALL conda env create -f ncfiles/vCapiPy.yml 
ECHO Environment creation finished!
ECHO ____________________________________________________________________________
TIMEOUT /t 3 /nobreak > NUL

ECHO Step 2/3: File modification to ensure CapiPy's correct operation...
MOVE %~dp0\ncfiles\CapiPy.bat %~dp0\ > NUL
SET /p LICENSE="Enter your MODELLER license (https://salilab.org/modeller/registration.html):"
CALL python ./ncfiles/file_modification.py %1 %LICENSE% > CON
ECHO File modification finished...
TIMEOUT /t 3 /nobreak > NUL
ECHO ____________________________________________________________________________
CALL tar.exe -cf CapiPy.zip CapiPy

ECHO Installation complete. Check both log files for errors.
pause