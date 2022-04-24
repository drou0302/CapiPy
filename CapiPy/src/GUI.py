import PySimpleGUI as sg
import sys
import os
import traceback
import src.BMOD as BMOD
import src.QUAT as QUAT
import src.AS_ID as AS_ID
import src.S2_analysis as S2_analysis
import CapiPy.ClusterDistance as ClusterDistance
import CapiPy.Ref_retrieval as Ref_retrieval
import webbrowser
import json
import subprocess
import time

from contextlib import redirect_stdout, redirect_stderr
from os import path
from sys import platform



def exes_location():
    if "linux" in platform:
        pymol_exe = "pymol"
        clustalw_exe = "clustalw"
        text_exe =  "gedit"
        csv_exe = 'soffice -calc'

    elif platform == "darwin":
        pymol_exe = r'/Applications/Pymol.app'
        clustalw_exe = r'/Applications/clustalw2'
        text_exe = r'/Applications/TextEdit.app'
        csv_exe = r'/Applications/Microsoft Excel.app'

    else:
        pymol_exe = r'C:/PyMOL/PyMOLWin.exe'
        clustalw_exe = r'C:/Program Files (x86)/ClustalW2/clustalw2.exe'
        text_exe = r'C:/Windows/system32/notepad.exe'
        csv_exe = r'C:/Program Files/Microsoft Office/root/Office16/EXCEL.EXE'

    return pymol_exe, clustalw_exe, text_exe, csv_exe

pymol_exe, clustalw_exe, text_exe, csv_exe = exes_location()
dir = BMOD.__file__.split('\\')
currentdir = ''
i=0
for p1 in dir:
    if i<len(dir)-2:
        currentdir += p1 + '/'
        i+=1

ex_loc = os.path.realpath(__file__).replace('GUI.py', '')
settings_file =  ex_loc + 'ncfiles/settings_file.cfg'
default_settings = {'Pymol_location': pymol_exe ,
                    'Clustalw_location': clustalw_exe,
                    'text_exe': text_exe,
                    'csv_exe': csv_exe,
                    'theme': 'Reddit'}
key_win = {'Pymol_location': '-PYMOL-',
           'Clustalw_location': '-CLUSTALW-',
           'text_exe': '-TEXT-',
           'csv_exe': '-CSV-',
           'theme': '-THEME-'}
settings = {}

if os.path.isfile(settings_file):
    settings = json.load(open(settings_file))
else:
    sg.popup_quick_message('No settings file found. Creating one with default values...',
                           keep_on_top=True, background_color='red', text_color='white')
    settings = default_settings
    with open(settings_file, "w+") as s:
        json.dump(default_settings, s)


def main():
    if 'win' in platform:
        scrpath = os.path.realpath(__file__).split('\\')
    else:
        scrpath = os.path.realpath(__file__).split('/')
    imgpath = ''
    for i in scrpath:
        if i != 'GUI.py':
            imgpath += str(i) + '/'
    img = imgpath + 'ncfiles/CapiPy.png'
    sg.theme(settings['theme'])

    layout_main = [[sg.Text('Welcome to CapiPy!', font =('Helvetica', 20, 'bold'))],
                   [sg.Button('1. Blast and Modeller', font =('Helvetica', 12), size=(30, 3))],
                   [sg.Button('1.1. Quaternary assembly', font =('Helvetica', 11), size=(30, 2))],
                   [sg.Button('2. Active site ID', font =('Helvetica', 12), size=(30, 3))],
                   [sg.Button('3. Surface and Clusters', font =('Helvetica', 12), size=(30, 3))],
                   [sg.Button('3.1. Cluster distance', font =('Helvetica', 11), size=(30, 2))],
                   [sg.Button('4. Immobilization papers retrieval', font =('Helvetica', 12), size=(30, 3))],
                   [sg.Image(img)],
                   [sg.Text('David Roura Padrosa, 2020', font =('Helvetica', 6))],
                   [sg.Button('Settings', size=(7, 1)), sg.Button('Help', size=(7, 1))]]
  
    window_main = sg.Window('CapiPy', layout_main, element_justification='c') 
    while True: 
        event_main, values_main = window_main.read() 
        if event_main in (None, 'Exit'): 
            break
        elif event_main == "1. Blast and Modeller":
            layout_m1 = [[sg.Text('Project name:'), sg.Input(key='-PROJECT-')],
                         [sg.Column([[sg.Text('Folder'), sg.In(size=(25,1), 
                         key='-FOLDER-'),sg.FolderBrowse()]], element_justification='c')],
                         [sg.Text('Please enter your protein sequence: ', font =('Helvetica', 10, 'bold'))],
                         [sg.Multiline(size=(100, 10), key='-SEQUENCE-')],
                         [sg.Text('How do you want to run BLAST?', font =('Helvetica', 10, 'bold'))],
                         [sg.Checkbox('LOCAL'), sg.Checkbox('WEB')],
                         [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')),
                          sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                         [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                         [sg.Output(key='-output-', size=(100, 20))],
                         [sg.Button('Run'), sg.Button('Exit')]]
            window_m1 = sg.Window('Blast and Modeller', layout_m1)
            while True: 
                event_m1, values_m1 = window_m1.read()
                checks = []
                if event_m1 in (None, 'Exit'):
                    break
                elif event_m1 == 'Run':
                    window_m1['-output-'].Update('')
                    
                    project_name = values_m1['-PROJECT-']
                    directory_name = values_m1['-FOLDER-'] 
                    folder_name = directory_name + '/' + project_name
                    
                    if os.path.exists(folder_name) == False:
                        os.mkdir(folder_name)
                    
                    if len(values_m1['-SEQUENCE-']) > 1:
                        prot_sequence = values_m1['-SEQUENCE-']
                    else:
                        checks.append("Error: enter your protein sequence!")
                    if values_m1[0] is True:
                        blast_response = "LOCAL"
                    elif values_m1[1] is True:
                        blast_response == "WEB"
                    else:
                        checks.append("Error: select how you want to run BLAST!")
                    if values_m1[2] is True:
                        output_programs = "YES"
                    elif values_m1[3] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")
                if len(checks) == 0:
                    print('Module running, be patient. This can take a couple of minutes...\n')
                    redirect_stdout(window_m1['-output-'])
                    window_m1.refresh()
                    try:
                        BMOD.main(folder_name, prot_sequence, blast_response, clustalw_exe)
                        if output_programs == "YES":
                            if os.path.isfile(folder_name + '/Blast&Modeller/query_multimeric.pdb'):
                                subprocess.Popen([settings['Pymol_location'],
                                                  folder_name + '/Blast&Modeller/query_monomer.pdb',
                                                  folder_name + '/Blast&Modeller/query_multimeric.pdb',
                                                  folder_name + '/Blast&Modeller/template.pdb'])
                            elif not os.path.isfile(folder_name + '/Blast&Modeller/query_multimeric.pdb'):
                                subprocess.Popen([settings['Pymol_location'], folder_name + '/Blast&Modeller/query_monomer.pdb',
                                                  folder_name + '/Blast&Modeller/template.pdb'])
                            else:
                                print("It seems CapiPy cannot find the template pdb file!")
                        os.chdir(directory_name)
                    except BaseException as error:
                        print('Something went wrong...\n' + str(error))
                        os.chdir(directory_name)
                else:
                    for sent in checks:
                        print(sent)
                redirect_stdout(window_m1['-output-'])
                window_m1.refresh()
            window_m1.close()    
       
        elif event_main == "1.1. Quaternary assembly":
            layout_m1_1 = [[sg.Text('Select your monomer pdb file:'), sg.FileBrowse(key='-FILE_MONOMER-')],
                           [sg.Text('Please indicate the template PDB code: ', 
                                    font =('Helvetica', 10, 'bold')), sg.Input(key='-PDB_CODE-')],
                           [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')),
                            sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                           [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                           [sg.Output(key='-output-', size=(80, 20))],
                           [sg.Button('Run'), sg.Button('Exit')]]
            window_m1_1 = sg.Window('Quaternary assembly', layout_m1_1)
            while True: 
                event_m1_1, values_m1_1 = window_m1_1.read()
                if event_m1_1 in (None, 'Exit'):
                    break
                elif event_m1_1 == 'Run':
                    window_m1_1['-output-'].Update('')
                    checks = []
                    if values_m1_1['-PDB_CODE-'] != '':
                        monomeric_model = values_m1_1['-FILE_MONOMER-']
                        folder_name = os.path.split(monomeric_model)[0]
                    else:
                        checks.append("Error: indicate the monomeric pdb!")
                    if values_m1_1['-PDB_CODE-'] != '':
                        multimeric_model = values_m1_1['-PDB_CODE-']
                    else:
                        checks.append("Error: indicate the template as a PDB code (four letter code)")
                    if values_m1_1[0] is True:
                        output_programs = "YES"
                    elif values_m1_1[1] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")
                    if len(checks) == 0:
                        print('Module running, be patient. This can take a couple of minutes...\n')
                        redirect_stdout(window_m1_1['-output-'])
                        window_m1_1.refresh()
                        try:
                            location_quat_pdb = QUAT.main(folder_name, monomeric_model, multimeric_model, clustalw_exe)
                            if output_programs == "YES":
                                print(location_quat_pdb)
                                if os.path.isfile(location_quat_pdb + '/quat_as_' +
                                              multimeric_model + '.pdb'):
                                    print("one")
                                    subprocess.Popen([settings['Pymol_location'],
                                                 location_quat_pdb + '/quat_as_' +
                                                  multimeric_model + '.pdb',
                                                 location_quat_pdb +'/quat_template.pdb'])
                                elif os.path.isfile(location_quat_pdb + '/Quaternary/quat_as_' +
                                    multimeric_model + '.pdb'):
                                    print("two")
                                    subprocess.Popen([settings['Pymol_location'],
                                                 location_quat_pdb + '/Quaternary/quat_as_' +
                                                  multimeric_model + '.pdb',
                                                 location_quat_pdb +'/Quaternary/quat_template.pdb'])
                                elif not os.path.isfile('/Blast&Modeller/quat_as_' +
                                    multimeric_model + '.pdb'):
                                    print("three")
                                    subprocess.Popen([settings['Pymol_location'],
                                                  location_quat_pdb + '/query_multimeric.pdb',
                                                  location_quat_pdb +'/quat_template.pdb'])
                                else:
                                    print("It seems CapiPy cannot find the pdb files!")
                            os.chdir(folder_name)
                        except BaseException as error:
                            print('Something went wrong...\n' + str(error))
                            os.chdir(folder_name)
                    else:
                        for sent in checks:
                            print(sent)
                    redirect_stdout(window_m1_1['-output-'])
                    window_m1_1.refresh()                

            window_m1_1.close()

        elif event_main == "2. Active site ID":
            layout_m2 = [[sg.Column([[sg.Text('Folder'), sg.In(size=(25,1), 
                         key='-FOLDER-'),sg.FolderBrowse()]], element_justification='c')],
                         [sg.Text('Have you run the Blast&Modeller for your query protein?',
                                  font =('Helvetica', 10, 'bold'))],
                         [sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                         [sg.Text('\tIf you have not run the Blast&Modeller or it is not in the same folder, \n\tplease indicate the protein sequence and the model protein to use:')],
                         [sg.Text('\tProtein sequence:'), sg.Multiline(size=(70, 10), key='-SEQUENCE-')],
                         [sg.Text('\tModel protein PDB code -4 characters identifier-:'), sg.Input(key='-MODEL-')],
                         [sg.Text('If available, do you want to use the data in UniProt for the estimation?',
                                  font =('Helvetica', 10, 'bold'))],
                         [sg.Checkbox('Yes, use UniProt preferentially.', default=True),
                          sg.Checkbox('No, only M-CSA data.')],
                         [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')),
                          sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                         [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                         [sg.Output(key='-output-', size=(95, 20))],
                         [sg.Button('Run'), sg.Button('Exit')]]
            window_m2 = sg.Window('Active Site ID', layout_m2)
            while True: 
                event_m2, values_m2 = window_m2.read()
                if event_m2 in (None, 'Exit'):
                    break
                elif event_m2 == 'Run':
                    window_m2['-output-'].Update('')
                    checks = []
                    folder_name = values_m2['-FOLDER-']
                    if values_m2[0] is True and values_m2[1] is False:
                        folder_answer = "NO"
                        prot_sequence = ''
                        model_prot = ''
                    elif values_m2[0] is False and values_m2[1] is True:
                        folder_answer = "YES"
                        if len(values_m2['-SEQUENCE-']) > 0:
                            prot_sequence = values_m2['-SEQUENCE-'].replace('\n', '')
                        if len(values_m2['-MODEL-']) == 4:
                            if values_m2['-MODEL-'].isalpha() is True:
                                checks.append("Error: It seems that your format is not a PDBid.")
                            elif values_m2['-MODEL-'].isnumeric() is True:
                                checks.append("Error: Remember that PDBid are formed by letters and numbers.")
                            else:
                                model_prot = values_m2['-MODEL-']
                        else:
                            checks.append('Error: Enter the protein sequence or select the query '
                                          'from a previously created folder!')
                    elif values_m2[0] == values_m2[2]:
                        checks.append("Error: create a folder or indicate its location!")

                    if values_m2[2] is True:
                        use_uniprot = 'YES'
                    elif values_m2[3] is True:
                        use_uniprot = 'NO'
                    else:
                        checks.append("Error: Select if you want to use UniProt preferentially or not.")

                    if values_m2[4] is True:
                        output_programs = "YES"
                    elif values_m2[5] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")

                    if len(checks) == 0:
                        print('Module running, be patient. This can take a couple of minutes...\n')
                        redirect_stdout(window_m2['-output-'])
                        window_m2.refresh()
                        try:
                            uniprotID = AS_ID.main(folder_name, folder_answer, prot_sequence,
                                                   model_prot, use_uniprot, clustalw_exe)
                            if output_programs == 'YES':
                                webbrowser.open('https://www.uniprot.org/uniprot/' + str(uniprotID))
                                time.sleep(3)
                                subprocess.Popen([settings['text_exe'],
                                                 folder_name + '/ActiveSite/Active_site.txt'])
                                time.sleep(3)
                                subprocess.Popen([settings['text_exe'],
                                                  folder_name + '/ActiveSite/refsequp.aln'])
                            os.chdir(folder_name)
                        except BaseException as error:
                            print('Something went wrong...\n' + str(error))
                            os.chdir(folder_name)
                    else:
                        for sent in checks:
                            print(sent)
                    redirect_stdout(window_m2['-output-'])
                    window_m2.refresh()                

            window_m2.close()

        elif event_main == "3. Surface and Clusters":
            layout_m3 = [[sg.Column([[sg.Text('Select the PDB file:'), sg.In(size=(25,1), 
                         key='-PDB_FILE-'),sg.FileBrowse()]], element_justification='c')],
                         [sg.Text('Active Site residues (Ex. 1, 2, 3):', font =('Helvetica', 10, 'bold'))],
                         [sg.Checkbox('ActiveSiteID results', default=True),
                          sg.Text('Or indicate residues:'), sg.Input()],
                         [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                         [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')),
                          sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                         [sg.Output(key='-output-', size=(80, 20))],
                         [sg.Button('Run'), sg.Button('Exit')]]
            window_m3 = sg.Window('Surface and Clusters', layout_m3)
            while True: 
                event_m3, values_m3 = window_m3.read()
                if event_m3 in (None, 'Exit'):
                    break
                elif event_m3 == 'Run':
                    window_m3['-output-'].Update('')
                    checks = []
                    pdbfile = values_m3['-PDB_FILE-'] 
                    folder_name = os.path.split(pdbfile)[0]
                    if values_m3[0] is True:
                        active_site = 'YES'
                    elif values_m3[1] is False and len(values_m3[6]) > 0:
                        active_site = values_m3[2].replace(" ", "").split(",")
                    else:
                        checks.append("Error: Untick the box or indicate the active site residues!")
                    if values_m3[2] is True:
                        output_programs = "YES"
                    elif values_m3[3] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")
                        
                    if len(checks) == 0:
                        print('Module running, be patient. This can take a few minutes...\n')
                        redirect_stdout(window_m3['-output-'])
                        window_m3.refresh()
                        try:
                            S2_analysis.main(folder_name, pdbfile, active_site)
                            if output_programs == 'YES':
                                try:
                                    subprocess.Popen([settings['Pymol_location'],
                                             folder_name + '/Size&Clusters/query.pdb',
                                             folder_name + '/Size&Clusters/PyMol_clusters.pml'])
                                    subprocess.Popen([settings['text_exe'],
                                             folder_name + '/Size&Clusters/General_info.txt'])
                                except BaseException:
                                    print("It seems CapiPy cannot find the template pdb file!")
                            os.chdir(folder_name)
                        except BaseException:
                            print('Something went wrong...\n' + traceback.format_exc())
                            os.chdir(folder_name)
                    else:
                        for sent in checks:
                            print(sent)
                    redirect_stdout(window_m3['-output-'])
                    window_m3.refresh() 

            window_m3.close()
       
        elif event_main == "3.1. Cluster distance":
         
            layout_m3_1 = [[sg.Column([[sg.Text('Select the PDB file:'), sg.In(size=(25,1), 
                          key='-PDB_FILE-'),sg.FileBrowse()]], element_justification='c')],
                          [sg.Text('Residues to calculate distance from (Ex. 1_A, 2_A, 300_B:', 
                                    font =('Helvetica', 10, 'bold'))],
                          [sg.Input()],
                          [sg.Text('Enter the file name you want to save the results (no spaces!):',
                          font =('Helvetica', 10, 'bold'))],
                          [sg.Input()],
                          [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')),
                           sg.Checkbox('Yes', default=True), sg.Checkbox('No')],
                          [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                          [sg.Output(key='-output-', size=(80, 20))],
                          [sg.Button('Run'), sg.Button('Exit')]]
        
            window_m3_1 = sg.Window('Cluster distance', layout_m3_1)
        
            while True: 
                event_m3_1, values_m3_1 = window_m3_1.read()
                if event_m3_1 in (None, 'Exit'):
                    break
                elif event_m3_1 == 'Run':
                    window_m3_1['-output-'].Update('')
                    checks = []
                    if len(values_m3_1[0]) > 0:
                        pdbfile = values_m3_1['-PDB_FILE-']
                        folder_name = os.path.split(pdbfile)[0]
                    else:
                        checks.append('Error: Please enter the name of the folder you want to process the data from.')
                    if len(values_m3_1[0]) > 0:
                        residues = values_m3_1[0].replace(" ", "").split(",")
                        for res in residues:
                            if res.isdigit() is True:
                                checks.append("\nError: Please, indicate the chain also!")
                            elif res.isalpha() is True:
                                checks.append("Error: There should be at least one number your residues.")
                            elif res == "" or res == " ":
                                checks.append("\nError: Please, enter at least one position.")
                    elif len(values_m3_1[0]) > 0:
                        checks.append('Error: Enter at least one residue identified with it\'s number and chain (Ex. 1_A)')
                
                    if len(values_m3_1[1]) > 0:
                        file_name = values_m3_1[1]
                    else:
                        checks.append('Error: Enter the file name, please')

                    if values_m3_1[2] is True:
                        output_programs = "YES"
                    elif values_m3_1[3] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")

                    if len(checks) == 0:
                        print('Module running, be patient. This can take a few minutes...\n')
                        redirect_stdout(window_m3_1['-output-'])
                        window_m3_1.refresh()
                        try:
                            pml_file = ClusterDistance.main(folder_name, pdbfile, residues, file_name)
                            print(pml_file)
                            if output_programs == "YES":
                                try:
                                    subprocess.Popen([settings['Pymol_location'],
                                                 folder_name + '/Size&Clusters/query.pdb',
                                                 folder_name + '/Size&Clusters/PyMol_clusters.pml',
                                                 folder_name + '/Size&Clusters/' + pml_file])
                                except BaseException:
                                    print("It seems CapiPy cannot find the template pdb file!")
                            os.chdir(folder_name)
                        except BaseException as error:
                            print('Something went wrong...\n' + traceback.format_exc())
                            os.chdir(folder_name)
                    else:
                        for sent in checks:
                            print(sent)
                    redirect_stdout(window_m3_1['-output-'])
                    window_m3_1.refresh() 

            window_m3_1.close()
    
        elif event_main == "4. Immobilization papers retrieval":

            layout_m4 = [[sg.Text('Project name:'), sg.FolderBrowse(key='-FOLDER-')],
                         [sg.Text('Do you want to use the same query as in the other modules?',
                                  font =('Helvetica', 10, 'bold'))],
                         [sg.Checkbox('Yes'), sg.Checkbox('No')],
                         [sg.Text('If not, please copy your sequence here:')],
                         [sg.Multiline(size=(80, 5))],
                         [sg.Text('Do you want to open the results directly?', font =('Helvetica', 10, 'bold')), sg.Checkbox('Yes', default=True),
                          sg.Checkbox('No')],
                         [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                         [sg.Output(key='-output-', size=(80, 20))],
                         [sg.Button('Run'), sg.Button('Exit')]]

            window_m4 = sg.Window('Imm. Papers retrieval', layout_m4)

            while True: 
                event_m4, values_m4 = window_m4.read()
                if event_m4 in (None, 'Exit'):
                    break
                elif event_m4 == 'Run':
                    window_m4['-output-'].Update('')
                    checks = []
                    folder_name = values_m4['-FOLDER-']
                    if values_m4[0] is True and values_m4[1] is False:
                        query = 'YES'
                    elif values_m4[0] is False and values_m4[1] is True:
                        if len(values_m4[3]) > 1:
                            query = values_m4[3].replace("\n", "")
                        else:
                            checks.append('Error: please enter a valid sequence. '
                                          'Remember, it has to be as one letter amino acid code!')
                    elif values_m4[0] == values_m4[1]:
                        checks.append('Error: Check ONE of the tickboxes!')

                    if values_m4[2] is True:
                        output_programs = "YES"
                    elif values_m4[3] is True:
                        output_programs = "NO"
                    else:
                        checks.append("Error: Select if you want your results to be shown directly or not.")

                    if len(checks) == 0:
                        print('Module running, be patient. This can take a few minutes...\n')
                        redirect_stdout(window_m4['-output-'])
                        window_m4.refresh()
                        try:
                            Ref_retrieval.main(folder_name, query)
                            if output_programs =="YES":
                                subprocess.Popen([settings['csv_exe'],
                                              folder_name + '/RelevantPapers/' + 'relevant_papers.csv'])
                            os.chdir(folder_name)
                        except BaseException:
                            print('Something went wrong...\n' + traceback.format_exc())
                            os.chdir(folder_name)
                    else:
                        for sent in checks:
                            print(sent)

                    redirect_stdout(window_m4['-output-'])
                    window_m4.refresh() 

            window_m4.close()
        
        elif event_main == 'Settings':
            sg.theme(settings['theme'])
            layout_set = [  [sg.Text('Settings', font='Any 12')],
                [sg.Text('Pymol location', font =('Helvetica', 10, 'bold')),
                 sg.Input(key='-PYMOL-'), sg.FileBrowse(target='-PYMOL-')],
                [sg.Text('Clustalw location', font =('Helvetica', 10, 'bold')),
                 sg.Input(key='-CLUSTALW-'), sg.FileBrowse(target='-CLUSTALW-')],
                [sg.Text('Text editor location', font =('Helvetica', 10, 'bold')),
                 sg.Input(key='-TEXT-'), sg.FileBrowse(target='-TEXT-')],
                [sg.Text('CSV reader location', font =('Helvetica', 10, 'bold')),
                 sg.Input(key='-CSV-'), sg.FileBrowse(target='-CSV-')],
                [sg.Text('Theme', font =('Helvetica', 10, 'bold')),
                 sg.Combo(sg.theme_list(), size=(20, 20), key='-THEME-')],
                [sg.Button('Save'), sg.Button('Exit')]]
            window_set = sg.Window('Settings', layout_set, finalize=True)
            
            for key in key_win:
                try:
                    window_set[key_win[key]].update(settings[key])
                except Exception:
                    print(f'Problem updating PySimpleGUI window from settings. Key = {key}')

            
            event_set, values_set = window_set.read()

            while True:
                if event_set in (None, 'Exit'):
                    break
                elif event_set == 'Save':
                    for key in key_win:  
                        try:
                            settings[key] = values_set[key_win[key]]
                        except Exception as e:
                            print(f'Error: Problem updating settings from window values. Key = {key}' + str(e))
                        with open(settings_file, 'w') as f:
                            json.dump(settings, f)
                    sg.popup('Settings saved')
                    break
            window_set.close()
        
        elif event_main == "Help":
            webbrowser.open('https://github.com/drou0302/CapiPy')
        else:
            continue
    window_main.close()

if __name__ == "__main__":
    main()
