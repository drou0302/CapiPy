import PySimpleGUI as sg
import sys
import os
import subprocess
from contextlib import redirect_stdout, redirect_stderr
import traceback
import BMOD, QUAT, AS_ID, S2_analysis,  ClusterDistance, Ref_retrieval
import webbrowser

location = os.getcwd()
sys.path.insert(1, location + '/Scripts/')

currentdir = os.getcwd()
img = currentdir + '/Installation/Environment/files/CapiPy.png'
sg.theme('Reddit') 
   
layout_main = [[sg.Text('Welcome to CapiPy!', font =('Helvetica', 20, 'bold'))],
               [sg.Button('1. Blast and Modeller', font =('Helvetica', 12), size=(30, 3))],
               [sg.Button('1.1. Quaternary assembly', font =('Helvetica', 12), size=(30, 3))],
               [sg.Button('2. Active site ID', font =('Helvetica', 12), size=(30, 3))],
               [sg.Button('3. Surface and Clusters', font =('Helvetica', 12), size=(30, 3))],
               [sg.Button('3.1. Cluster distance', font =('Helvetica', 12), size=(30, 3))],
               [sg.Button('4. Immobilization papers retreival', font =('Helvetica', 12), size=(30, 3))],
               [sg.Image(img)],
               [sg.Text('David Roura Padrosa, 2020', font =('Helvetica', 6))],
               [sg.Button('Help', size=(7, 1))]]
  
window_main = sg.Window('CapiPy', layout_main, element_justification='c') 
  
while True: 
    event_main, values_main = window_main.read() 
    if event_main in (None, 'Exit'): 
        break
    elif event_main == "1. Blast and Modeller":
        layout_m1 = [[sg.Text('Do you need to create a new folder?', font =('Helvetica', 10, 'bold'))], 
                     [sg.Checkbox('Yes'), sg.Text('Enter the name:'), sg.Input()],
                     [sg.Checkbox('No'), sg.Text('Enter the location:'), sg.Input()],
                     [sg.Text('Please enter your protein sequence: ', font =('Helvetica', 10, 'bold'))],
                     [sg.Multiline(size=(100, 10))],
                     [sg.Text('How do you want to run BLAST?', font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('LOCAL'), sg.Checkbox('WEB')],
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
                if values_m1[0] is True and values_m1[2] is False:
                    folder_answer = "YES"
                    folder_name = values_m1[1]
                elif values_m1[0] is False and values_m1[2] is True:
                    folder_answer = "NO"
                    folder_name = values_m1[3]
                elif values_m1[0] is False and values_m1[2] is False:
                    checks.append("Missing information: create a folder or indicate its location!")
                else:
                    checks.append("Missing information: create a folder or indicate its location, but not both!")
                if len(values_m1[4]) > 1:
                    prot_sequence = values_m1[4]
                else:
                    checks.append("Missing information: enter your protein sequence!")
                if values_m1[5] is True:
                    blast_response = "LOCAL"
                elif values_m1[6] is True:
                    blast_response == "WEB"
                else:
                    checks.append("Missing information: select how you want to run BLAST!")
            if len(checks) == 0:
                print('Module running, be patient. This can take a couple of minutes...\n')
                redirect_stdout(window_m1['-output-'])
                window_m1.refresh()
                try:
                    BMOD.main(folder_answer, folder_name, prot_sequence, blast_response)
                except BaseException as error:
                    print('Something went wrong...\n' + str(error))
            else:
                for sent in checks:
                    print(sent)
            redirect_stdout(window_m1['-output-'])
            window_m1.refresh()
        window_m1.close()    
       
    elif event_main == "1.1. Quaternary assembly":
        layout_m1_1 = [[sg.Text('Enter the location of your pdb file (monomer):', 
                                font =('Helvetica', 10, 'bold')), sg.Input()],
                       [sg.Text('Please indicate the monomeric model: ', 
                                font =('Helvetica', 10, 'bold')), 
                        sg.Input(default_text='query_monomer.pdb')],
                       [sg.Text('Please indicate the template PDB code: ', 
                                font =('Helvetica', 10, 'bold')), sg.Input()],
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
                if values_m1_1[0] != '':
                    folder_name = values_m1_1[0]
                else:
                    checks.append("Missing information: indicate the folder location!")
                if values_m1_1[1] != '':
                    monomeric_model = values_m1_1[1]
                else:
                    checks.append("Missing information: indicate the monomeric pdb!")
                if values_m1_1[2] != '':
                    multimeric_model = values_m1_1[2]
                else:
                    checks.append("Missing information: indicate the template as a PDB code (four letter code)")
                if len(checks) == 0:
                    print('Module running, be patient. This can take a couple of minutes...\n')
                    redirect_stdout(window_m1_1['-output-'])
                    window_m1_1.refresh()
                    try:
                        QUAT.main(folder_name, monomeric_model, multimeric_model)
                    except BaseException as error:
                        print('Something went wrong...\n' + str(error))
                else:
                    for sent in checks:
                        print(sent)
                redirect_stdout(window_m1_1['-output-'])
                window_m1_1.refresh()                

        window_m1_1.close()

    elif event_main == "2. Active site ID":
        layout_m2 = [[sg.Text('Do you need to create a new folder?', font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('No', default=True), sg.Text('Enter the folder location:'), sg.Input()],
                     [sg.Checkbox('Yes'), sg.Text('Enter the name:'), sg.Input()],
                     [sg.Text('\tIf you are creating a new folder, please indicate the protein sequence '
                              '- in one letter code - :')],
                     [sg.Text('\tProtein sequence:'), sg.Multiline(size=(70, 10))],
                     [sg.Text('\tModel protein PDB code -4 characters identifier-:'), sg.Input()],
                     [sg.Text('If available, do you want to use the data in UniProt for the estimation?',
                              font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('Yes, use UniProt preferentially.', default=True),
                      sg.Checkbox('No, only M-CSA data.')],
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
                if values_m2[0] is True and values_m2[2] is False:
                    folder_answer = "NO"
                    folder_name = values_m2[1]
                    prot_sequence = ''
                    model_prot = ''
                elif values_m2[0] is False and values_m2[2] is True:
                    folder_answer = "YES"
                    folder_name = values_m2[3]
                    if len(values_m2[4]) > 0:
                        prot_sequence = values_m2[4].replace('\n', '')
                    if len(values_m2[5]) == 4:
                        if values_m2[5].isalpha() is True:
                            checks.append("Error: It seems that your format is not a PDBid.")
                        elif values_m2[5].isnumeric() is True:
                            checks.append("Error: Remember that PDBid are formed by letters and numbers.")
                        else:
                            model_prot = values_m2[5]
                    else:
                        checks.append('Missing information: Enter the protein sequence or select the query '
                                      'from a previously created folder!')
                elif values_m2[0] == values_m2[2]:
                    checks.append("Missing information: create a folder or indicate its location!")

                if values_m2[5] is True:
                    use_uniprot = 'YES'
                elif values_m2[6] is True:
                    use_uniprot = 'NO'

                if len(checks) == 0:
                    print('Module running, be patient. This can take a couple of minutes...\n')
                    redirect_stdout(window_m2['-output-'])
                    window_m2.refresh()
                    try:
                        AS_ID.main(folder_name, folder_answer, prot_sequence, model_prot, use_uniprot)
                    except BaseException as error:
                        print('Something went wrong...\n' + str(error))
                else:
                    for sent in checks:
                        print(sent)
                redirect_stdout(window_m2['-output-'])
                window_m2.refresh()                

        window_m2.close()

    elif event_main == "3. Surface and Clusters":
        layout_m3 = [[sg.Text('Do you need to create a new folder?', font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('No', default=True), sg.Text('\tEnter the folder name:'), sg.Input()],
                     [sg.Checkbox('Yes'), sg.Text('\tGive a name to the folder:'), sg.Input()],
                     [sg.Text('Indicate the location of the pdb file to use:', 
                              font =('Helvetica', 10, 'bold'))],
                     [sg.Input(default_text='query_multimeric.pdb')],
                     [sg.Text('Active Site residues (Ex. 1, 2, 3):', font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('ActiveSiteID results', default=True), sg.Text('Or indicate residues:'), sg.Input()],
                     [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
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
                if values_m3[0] is True and values_m3[2] is False:
                    folder_answer = "NO"
                    if values_m3[1] != "":
                        folder_name = values_m3[1]
                    else:
                        checks.append("Error: Indicate the name of the folder.")
                        
                elif values_m3[0] is False and values_m3[2] is True:
                    folder_anwer = 'YES'
                    if values_m3[1] != "":
                        folder_name = values_m3[3]
                    else:
                        checks.append("Error: Indicate the name of the folder.")
                elif values_m3[0] == values_m3[2]:
                    checks.append('Error: Select a folder or create a new one.')   
                         
                if values_m3[4] != '':
                    pdbfile = values_m3[4]
                else:
                        checks.append('Error: Indicate the name of the pdb file!')
    
                if values_m3[5] is True:
                    active_site = 'YES'
                elif values_m3[5] is False and len(values_m3[6] > 0):
                    active_site = values_m3[6].replace(" ", "").split(",")
                else:
                    checks.append("Error: Untick the box or indicate the active site residues!")

                if len(checks) == 0:
                    print('Module running, be patient. This can take a few minutes...\n')
                    redirect_stdout(window_m3['-output-'])
                    window_m3.refresh()
                    try:
                        S2_analysis.main(folder_answer, folder_name, pdbfile, active_site)
                    except BaseException as error:
                        print('Something went wrong...\n' + traceback.format_exc())
                else:
                    for sent in checks:
                        print(sent)
                redirect_stdout(window_m3['-output-'])
                window_m3.refresh() 

        window_m3.close()
       
    elif event_main == "3.1. Cluster distance":
         
        layout_m3_1 = [[sg.Text('Enter the folder name:', font =('Helvetica', 10, 'bold'))],
                      [sg.Input()],
                      [sg.Text('Residues to calculate distance from (Ex. 1_A, 2_A, 300_B:', 
                                font =('Helvetica', 10, 'bold'))],
                       [sg.Input()],
                       [sg.Text('Enter the file name you want to save the results (no spaces!):',
                        font =('Helvetica', 10, 'bold'))], 
                       [sg.Input()],
                       [sg.Text('OUTPUT:', font =('Helvetica', 10, 'bold'))],
                       [sg.Output(key='-output-', size=(80, 20))],
                       [sg.Button('Run'), sg.Button('Exit')]]
        
        window_m3_1 = sg.Window('Cluster distance', layout_m3_1)
        
        while True: 
            event_m3_1, values_m3_1 = window_m3_1.read()
            print(event_m3_1, values_m3_1)
            if event_m3_1 in (None, 'Exit'):
                break
            elif event_m3_1 == 'Run':
                window_m3_1['-output-'].Update('')
                checks = []
                if len(values_m3_1[0]) > 0:
                    folder_name = values_m3_1[0]
                else:
                    checks.append('Error: Please enter the name of the folder you want to process the data from.')
                if len(values_m3_1[1]) > 0:
                    residues = values_m3_1[1].replace(" ", "").split(",")
                    for res in residues:
                        if res.isdigit() is True:
                            checks.append("\nPlease, indicate the chain also!")
                        elif res.isalpha() is True:
                            checks.append("There should be at least one number in there.")
                        elif res == "" or res == " ":
                            checks.append("\nPlease, enter at least one position.")
                elif len(values_m3_1[1]) > 0:
                    checks.append('Error: Enter at least one residue identified with it\'s number and chain (Ex. 1_A)')
                
                if len(values_m3_1[2]) > 0:
                    file_name = values_m3_1[2]
                else:
                    checks.append('Error: Enter the file name, please')
                if len(checks) == 0:
                    print('Module running, be patient. This can take a few minutes...\n')
                    redirect_stdout(window_m3_1['-output-'])
                    window_m3_1.refresh()
                    try:
                        ClusterDistance.main(folder_name, residues, file_name)
                    except BaseException as error:
                        print('Something went wrong...\n' + traceback.format_exc())
                else:
                    for sent in checks:
                        print(sent)
                redirect_stdout(window_m3_1['-output-'])
                window_m3_1.refresh() 

        window_m3_1.close()
    
    elif event_main == "4. Immobilization papers retreival":
        layout_m4 = [[sg.Text('Enter the folder name. You can enter an already existing folder name or '
                              'create a new one.', font =('Helvetica', 10, 'bold'))],
                     [sg.Input()],
                     [sg.Text('Do you want to use the same query as in the other modules?',
                              font =('Helvetica', 10, 'bold'))],
                     [sg.Checkbox('Yes'), sg.Checkbox('No')],
                     [sg.Text('If not, please copy your sequence here:')],
                     [sg.Multiline(size=(80, 5))],
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
                if len(values_m4[0]) > 0:
                    folder_name = values_m4[0]
                else:
                    checks.append('Error: Please enter the name of the folder you want to process the data from.')

                if values_m4[1] is True and values_m4[2] is False:
                    query = 'YES'
                elif values_m4[1] is False and values_m4[2] is True:
                    if len(values_m4[3]) > 1:
                        query = values_m4[3].replace("\n", "")
                    else:
                        checks.append('Error: please enter a valid sequence. '
                                      'Remember, it has to be as one letter amino acid code!')
                elif values_m4[1] == values_m4[2]:
                    checks.append('Error: Check ONE of the tickboxes!')

                if len(checks) == 0:
                    print('Module running, be patient. This can take a few minutes...\n')
                    redirect_stdout(window_m4['-output-'])
                    window_m4.refresh()
                    try:
                        Ref_retrieval.main(folder_name, query)
                    except BaseException as error:
                        print('Something went wrong...\n' + traceback.format_exc())
                else:
                    for sent in checks:
                        print(sent)

                redirect_stdout(window_m4['-output-'])
                window_m4.refresh() 
        window_m4.close()     
    elif event_main == "Help":
        webbrowser.open('https://github.com/drou0302/CapiPy')
    else:
        continue

window_main.close()

layout = [[sg.Text('Welcome to CapiPy-GUI!')],
          [sg.Text('Please enter your protein sequence: ')], 
          [sg.Input(key='-IN-')], 
          [sg.Button('Run'), sg.Button('Exit')]]
