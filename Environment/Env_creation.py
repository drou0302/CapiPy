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

print("Starting creation of vCapiPy environment...")

sys.stdout = open("install.log", "w")

cmd = "conda env create -f \"" + path + "/vCapiPy.yml\""
os.system(cmd)
sys.stdout.close()

print("All done, check the install.log file for errors!")
