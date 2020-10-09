#!/bin/bash
echo "Removing vCapiPy environment..."
conda env remove --name vCapiPy --all > uninstall.log 2>&1
rm ../*.sh
echo "Uninstall completed. Check uninstall.log for more information."
