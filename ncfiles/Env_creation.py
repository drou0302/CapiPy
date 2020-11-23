#! python3

# Installation directly with a python script"
import os
import subprocess
import time
import sys

# Create environment with necessary packages
full_path = os.path.realpath(__file__)
print("Locating the .yml file...")
path, filename = os.path.split(full_path)

log = open("ENV-CREATION.log", "w+")
print("Starting creation of vCapiPy environment...")
print("Starting creation of vCapiPy environment...", file="ENV-CREATION.log")

cmd = "conda env create -f \"" + path + "/vCapiPy.yml\""
sys.stdout, sys.stderr = log
os.system(cmd)

print("All done, check the install.log file for errors!")
