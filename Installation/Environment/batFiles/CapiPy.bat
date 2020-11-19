@echo off
call conda activate vCapiPy
python ./CapiPy/GUI.py
call conda deactivate
