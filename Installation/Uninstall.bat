@echo off
echo "Starting with environment deletion...""
call conda env remove --name vCapiPy >> Uninstall.log 2>&1
del ../*.bat 
echo "Uninstall complete. Check Uninstall.log for more information."
pause