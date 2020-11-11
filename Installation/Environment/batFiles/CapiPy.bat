@echo off
call conda activate vCapiPy
python ./Scripts/CapiPy.py
call conda deactivate
