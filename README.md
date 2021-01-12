# CapiPy
## Table of contents:
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Installation](#installation)
    * [Configuration](#configuration)
* [Troubleshooting](#troubleshooting)

## Introduction: 
__What is CapiPy?__
CapiPy (Computer Assistance for Protein Immobilisation – Python) is a small collection of 4 main packages to retrieve useful information for the immobilisation of a query protein.	 
Information from each module, specially the Blast and modeller is preferentially used by the others. In case no such information is given, the user will have to input the necessary information.

__MODULE 1: BLAST and MODELLER:__ this module only requires the 1 letter amino acid sequence of the query protein. This will be blasted against the protein data base to identify the protein with the highest similarity and a model will be created. If the protein is multimeric, the option to model the quaternary structure is also given via superposition of the monomeric model onto the other chains.

__MODULE 2: ACTIVE SITE ID:__ the identification of the active site from the query protein is performed using the information present in UNIPROT and/or the M-CSA database. If the protein is annotated to use any cofactor, this information is also given.

__MODULE 3: SIZE AND CLUSTERS IDENTIFICATION:__ this is the most important module regarding immobilisation. First, the surface and volume of the model of the query protein (or a user specified PDB) are calculated. Secondly, all protein residues are classified according to their solvent exposure and from the exposed ones, clusters formed by residues important for immobilisation identified. The output can easily be visualised in PyMOL.

__MODULE 4: IMMOBILISATION PAPER RETRIEVAL:__ to complement the information given by the other modules, in this last module 20 immobilisation papers related to the query protein are retrieved. To do so, the query used in the previous modules (or input here by the user) is blasted against the SwissProt database to identify the best hit.  The information of this best hit, such as protein name or identifiers, is extracted and used as keywords to submit a search into the Pubmed database. The 20 most relevant publications (by citation) are listed in an excel file.  

In addition, two stand-alone functionalities of the first and third module are also available:

__MODULE 1.1:__ QUATERNARY STRUCTURE DETERMINATION: with a monomeric model created, this module allows the user to select a template from the PDB dataset to create a different quaternary assembly.

__MODULE 3.1:__ CLUSTER DISTANCE: once the clusters have been identified, similar to the last part of the third module, this allows the calculation of the distance between the clusters and any user-specified position in the query protein.

## Technologies:
External software requirements:
- Anaconda (miniconda or anaconda) - https://www.anaconda.com/products/individual
- PyMOL - https://pymol.org/2/ or https://sourceforge.net/projects/pymol/files/latest/download
- Visual studio C++ tools - https://visualstudio.microsoft.com/downloads/
- Clustalw - http://clustal.org/clustal2/
- Blast 2.XX - https://www.ncbi.nlm.nih.gov/books/NBK279671/
- Modeller - https://salilab.org/modeller/
Included in the standard installation:
- Python 3.7 or later including:
    - biopython 1.77
    - metapub 0.5.5
    - more-itertools 8.4.0
    - numpy 1.19.0
    - PySimpleGUI 4.29.0
## Installation
- __Option 1: Create a separate environment using Anaconda - EASY INSTALLATION__
    Donwload the content of this repository and unpack the downloaded zip file.
    * If your OS is Windows, double click on the Install_CapiPy.bat and follow the instructions.
    * If your OS is Linux-based or MacOS, open a terminal window and type:
    ```
    $ cd ../CapiPy-main
    $ chmod +x Install_CapiPy.sh
    $ ./Install_CapiPy.sh
    ```
    _Note: During the installation, a zip file containing the source code will be created. If you modify any of the code in the CapiPy folder, remember to update the zip file also. Otherwise, the changes would not have any effect._
    ### Running CapiPy:
    To run CapiPy:
    - If your OS is Windows, double click on the CapiPy.bat
    - If your OS is Linux-based or MacOS:
    ```
    $ cd ../CapiPy-main
    $ chmod u+x CapiPy.sh
    $ ./CapiPy.sh
    ```
    or double click on the CapiPy.sh file and select open with Terminal.

- __Option 2: Add CapiPy as a package in your current python installation:__
    
    If you don't want to create a new environment, which is recomended, you can install CapiPy along with the necessary packages using pip. _Make sure you have            python 3.6 or later!_
    - To check you python version, open a terminal and type:
    ```
    $ python -V
    ```
    If your version is python 3.6 or later, continue with the installation. If not, please update your python. 
    - Once you have Python 3.6 or later, run the following command to install pip (if not installed already):
    ```
    $ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    $ python get-pip.py
    ```
    - To install the necessary packages, type:
    ```
    $ pip install biopython metapub more-itertools PySimpleGUI
    ```
    - And finally, to install CapiPy, enter:
    ```
    $ pip install CapiPy
    ```
     ### Running CapiPy:
    To run CapiPy:
    ```
    $ python -m CapiPy
    ```
### Configuration 

On the first use, make sure to check the configuration tab to tell CapiPy where to find PyMOL, ClustalW, a generic text editor and a CSV reader.
The default values are:
- Windows: 
    - C:/PyMOL/PyMOLWin.exe
    - C:/Program Files (x86)/ClustalW2/clustalw2.exe
    - C/Windows/system32/notepad.exe
    - C/Windows/system32/notepad.exe
    - C:/Program Files/Microsoft Office/root/Office16/EXCEL.EXE
- MacOS:
    - /Applications/Pymol.app
    - /Applications/clustalw2
    - /Applications/TextEdit.app
    - /Applications/Microsoft Excel.app
- Linux:
    - pymol
    - clustalw
    - gedit
    - soffice -calc

## Troubleshooting

### Installation

| Error | Troubleshoot |
| ---|---|
| CondaValueError: Value error: prefix already exists: | An environment with the same name already exists. Please delete it either by running the Uninstall.bat / Uninstall.sh file or run in a terminal: ````conda env remove --name vCapiPy````|
| Installation fails with error related to Python Levenshtein   | Make sure you have installed Microsoft Visual Studio and downloaded the default C++ modules, including: MSVC build tools and Windows 10 SDK. If you are still having the error, install Python levenshtein from the wheel file directly. You can find those in https://www.lfd.uci.edu/~gohlke/pythonlibs/#python-levenshtein. If you are trying to install CapiPy in a conda environment open the terminal and type: ``` $ conda activate vCapiPy ; $ cd PATH/TO/python-Levenshtein-0.12.0-cpXX-cpXX-win_amd64.whl ; $ pip install python-Levenshtein-0.12.0-cpXX-cpXX-win_amd64.whl```|
|Cannot fix the Residue.py file|Copy manually the Residue.py provided in the Environment folder to PATHtoCONDA\envs\vCapiPy\lib\site-packages\Bio\PDB\|
|Cannot fix the \_\_init\_\_.py Modeller|Go to the modeller folder (PATHtoCONDA\envs\vCapiPy\Library\modeller\modlib\modeller\) and manually edit the \_\_init\_\_.py file in line 68 to read: dpath = config.install_dir + '\\modlib\\'  __Make sure you do not change the indentation!__. If that doesn't work or the \_\_init\_\_.py file doesn't have such information, make sure that the config.py file in the same directory in the first line points to the correct directory of modeller. |
| Cannot fix the config file from Modeller | Go to the modeller folder (PATHtoCONDA\envs\vCapiPy\Library\modeller\modlib\modeller\) and manually edit the config file. Replace the XXXX after license to your modeller license code. |
|ModuleNotFoundError: No module named 'modeller'| Make sure the modeller locations are part of the PYTHONPATH variable. (Check https://salilab.org/modeller/9.25/release.html)| 



### Running:
#### External software:
| Error | Troubleshoot |
| ---|---|
| blastp not recognized as an internal or external command.   | Check that the folder containing the executable (PATH/blast-2.XX-/bin) is in your PATH variables.  |
| Local BLAST search does not work | Check that the databases are in the correct location and the ncbi.ini file is in the blast/bin folder and it points to the database folder. More information here https://www.ncbi.nlm.nih.gov/books/NBK279695/|
| Web BLAST search does not work or takes long time | Try again in a few minutes. Online BLAST depends on the server availability. |
| Modeller license is missing | Check the config.py file in your modeller installation. In this file, the second line should read: license = r'XXXX', where XXXX should be your license.|
| Modeller cannot find the specified location | - Make sure your config.py file of your modeller installation, points out to the correct modeller directory (p.e.: install_dir = r'C:/Users/User/conda/envs/vCapiPy/Library/modeller' for a virtual environment in Windows.- |
| Clustalw2 not recognized | Make sure clustalw2 is in the Applications folder and this in the PATH variables.|


#### CapiPy related errors:
| Error | Troubleshoot |
| ---|---|
|Cannot execute .sh files | You need to make the .sh file executable first with the command ``` $ chmod +x scriptname.sh``` |
|Cannot create a file when that file already exists | You are trying to run CapiPy and create a folder with a name that already exists in your current working directory |
