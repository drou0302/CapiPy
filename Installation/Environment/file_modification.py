import os
from sys import platform
import sys
import shutil

if "linux" in platform or "darwin" in platform:
    searchpath = "."
else:
    searchpath = "C:\\"

# Search for vCapiPy location#
for root, dir, files in os.walk(searchpath):
    if "vCapiPy" in dir:
        environmentloc = os.path.join(root, "vCapiPy")

# Replace Bio.PDB Residue.py file#
for root, dir, files in os.walk(environmentloc):
    if "Residue.py" in files:
        Resreplace = os.path.join(root, "Residue.py")
shutil.copyfile("./Environment/Files/Residue.py", Resreplace)

print("Fixed Bio.PDB module.")

# Search for __init__.py file from Modeller
locations = []
for root, dir, files in os.walk(environmentloc):
    if "__init__.py" in files:
        locations.append(os.path.join(root, "__init__.py"))
for locs in locations:
    if "vCapiPy\\Library\\modeller\\modlib\\modeller\\__init__.py" in locs:
        modreplace = locs
fixedfile = modreplace + "1"
shutil.copyfile(modreplace, fixedfile)
with open(fixedfile, "w") as toedit:
    with open(modreplace) as f:
        for line in f:
            toedit.write(line.replace("dpath = config.install_dir + '\\lib\\%s' "
                                      "% exetype", "dpath = config.install_dir + '\\modlib\\'"))

os.remove(modreplace)
shutil.copyfile(fixedfile, modreplace)
os.remove(fixedfile)

# Search for config.py file from Modeller
license = sys.argv[1]
print(license)
locations1 = []
for root, dir, files in os.walk(environmentloc):
    if "config.py" in files:
        locations1.append(os.path.join(root, "config.py"))
for locs in locations1:
    if "vCapiPy\\Library\\modeller\\modlib\\modeller\\config.py" in locs:
        lic_file = locs
lic_replace = lic_file + "1"
shutil.copyfile(lic_file, lic_replace)
with open(lic_replace, "w") as toedit:
    with open(lic_file) as f:
        for line in f:
            toedit.write(line.replace("license = r'XXXX' ", "license = r'" + str(license) + "'"))
os.remove(lic_file)
shutil.copyfile(lic_replace, lic_file)
os.remove(lic_replace)

print("Fixed modeller init and config file.")