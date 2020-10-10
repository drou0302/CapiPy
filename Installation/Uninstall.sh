#!/bin/bash
echo "Indicate the full location of conda.sh"
read $CONDA 
source $CONDA
echo "Removing vCapiPy environment..."
conda env remove --name vCapiPy --all > uninstall.log 2>&1
rm ../*.sh
echo "Uninstall completed. Check uninstall.log for more information."
read -rsp $'Press enter to close the finish.\n'
