@echo off
ECHO "Starting with environment deletion..."
CALL conda env remove --name vCapiPy >> Uninstall.log 2>&1
MOVE CapiPy.bat %~dp0\ncfiles\
DEL *.log

echo "Uninstall complete. Check Uninstall.log for more information."
pause