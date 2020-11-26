#!/bin/bash
source activate vCapiPy
DIR=$(dirname "$0")
cd "$DIR"
python CapiPy.zip
conda deactivate
