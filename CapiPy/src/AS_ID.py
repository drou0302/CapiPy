#!/usr/bin/env python3

# CatSiteID v.0.1
# Copyleft David Roura, 2020
"""Active site identification by using the information in PDB, UniProt and M-CSA in combination

"""

import urllib.request
import urllib.parse
import re
import csv
import os
import time
import shutil
import glob
import warnings

from sys import platform
from shutil import copyfile
from Bio import Align, SeqIO, SearchIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import BiopythonWarning
from Bio.Align import substitution_matrices

warnings.simplefilter('ignore', BiopythonWarning)


def main(folder_name, folder_answer, prot_sequence, model_prot, use_uniprot, clustalw_exe):
    initial_location = os.getcwd()
    working_directory = folder_name
    if os.path.exists(working_directory) is True:
        os.chdir(working_directory)
        # Check if ActiveSite directory exists
        if os.path.exists("ActiveSite") is False:
            os.mkdir("ActiveSite")
            os.chdir("./ActiveSite")
    else:
        print("The folder specified doesn't exist! Check again.")
        quit()

    # Download M-CSA in csv format if no in there already - https://www.ebi.ac.uk/thornton-srv/m-csa/

    url_atlas = "https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv"
    print("Downloading the curated data from the Mechanism and Catalytic Site atlas - www.ebi.ac.uk -.")
    with urllib.request.urlopen(url_atlas) as data, open('./curated_data.csv', 'w', encoding="utf-8") as f:
        f.write(data.read().decode())

    # Define the query sequence. Ask to use either the one in the blast&modeller directory or specify the one to use.
    id_1 = ''
    if folder_answer == 'NO':
        if os.path.exists("../Blast&Modeller") is True:
            shutil.copy("../Blast&Modeller/query.fasta", "./query.fasta")
            shutil.copy("../Blast&Modeller/model.fasta", "./model.fasta")
            with open("./model.fasta") as file:
                id_1 = file.readline().replace(">", "").rstrip().upper()
    elif folder_answer == "YES":
        with open("query.fasta", "w+") as f:
            protein_sequence = prot_sequence
            f.write(">query" + "\n" + protein_sequence)
        id_1 = model_prot
    # From pdb id to uniprot id, used in the M-CSA file#
    print("Searching for the uniprot equivalent of your pdb file...")
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from': 'PDB_ID', 'to': 'ACC', 'format': 'tab', 'query': id_1}
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    try:
        with urllib.request.urlopen(req) as f:
            response = f.read()
    except BaseException:
        print("There has been some problem connecting to UniPROT. Check your internet connection and try again!")
        time.sleep(3)
        quit()

    fromtoid = response.decode('utf-8')
    try:
        uniprotfullID = fromtoid.split("\t")[2]
        uniprotID = uniprotfullID.split("\n")[0]
    except BaseException:
        print("It seems that your model protein does not have a uniprot equivalent... "
              "You can try with an alternate model!")
        time.sleep(3)
        quit()
    if uniprotID.isalpha():
        print("Your pdb does not appear to have a UniProt equivalent. Sorry!")
    else:
        print("Uniprot id found. Extracting information...")
    handle = urllib.request.urlopen("https://www.uniprot.org/uniprot/"+str(uniprotID)+".xml")
    record = SeqIO.read(handle, "uniprot-xml")

    # Extract basic info from the uniprot id from the pdb file#

    protein_type = ""
    proteinECnumber = ""

    for information in record:
        if "submittedName_fullName" in record.annotations:
            protein_type = record.annotations["submittedName_fullName"]
            break
        elif "recommendedName_fullName" in record.annotations:
            protein_type = str(record.annotations["recommendedName_fullName"])
            break

    for information in record:
        if "submittedName_ecNumber" in record.annotations:
            proteinECnumber = record.annotations['submittedName_ecNumber']
            break
        elif "recommendedName_ecNumber" in record.annotations:
            proteinECnumber = record.annotations['recommendedName_ecNumber']
            break

    if isinstance(proteinECnumber, list):
        proteinECnumber = str(proteinECnumber[0])

    if proteinECnumber == "":
        print("It seems the information provided by uniprot is not enough. "
              "Please, check https://www.uniprot.org/uniprot/"+str(uniprotID))
        print("Your protein doesn't seem to have an EC number... Not much we can do without it!")
        if proteinECnumber.upper().replace(" ", "") == "NO":
            print("Sorry, but without the EC number there is little we can do...")
            time.sleep(3)
            quit()

    # If uniprot provided enough data, print into the console the answers
    res = open("Active_site.txt", "a")
    if isinstance(protein_type, list) is True:
        print("\nYour model protein, " + str(id_1)+", has been identified as a " + protein_type[0] +
              " with EC number " + str(proteinECnumber) + ".\n")
        print("\nYour model protein, " + str(id_1)+", has been identified as a " + protein_type[0] +
              " with EC number " + str(proteinECnumber) + ".\n", file=res)
    elif isinstance(protein_type, str) is True:
        print("\nYour model protein, " + str(id_1)+", has been identified as a "
              + protein_type.replace("[", "").replace("]", "") + " with EC number " + str(proteinECnumber) + ".\n")
        print("\nYour model protein, " + str(id_1)+", has been identified as a "
              + protein_type.replace("[", "").replace("]", "") + " with EC number " + str(proteinECnumber) + ".\n",
              file=res)
    else:
        print("\nYour model protein, " + str(id_1)+", has EC number " + str(proteinECnumber) + ".\n")
        print("\nYour model protein, " + str(id_1)+", has EC number " + str(proteinECnumber) + ".\n", file=res)
    res.close()
    time.sleep(3)

    # In a few cases, uniprot provides already the information about the catalytic site. In this case, we will compare
    # the pdb file provided to the query and avoid using the ATLAS site.
    act_site = []
    bind_site = []
    catalytic_res = {}
    cofactor = []

    # Find if information about the catalytic residues exists in the uniprot xml
    for i in range(len(record.features)):
        if record.features[i].type == "active site":
            act_site.append(int(record.features[i].location.end)-1)
        elif record.features[i].type == "binding site":
            bind_site.append(int(record.features[i].location.end)-1)

    cofactorslist = ["NAD", "FAD", "FMN", "SAM", "PLP", "ATP", "UTP", "ADP", "UDP", "CTP", "PAPS", "acetyl COA", "Zn",
                     "Fe", "Cu", "K", "Mg", "Mo", "Ni", "Se"]
    if "comment_cofactor" in record.annotations:
        for ncof in cofactorslist:
            for act in range(len(record.annotations["comment_cofactor"])):
                if ncof in record.annotations["comment_cofactor"][act]:
                    cofactor.append(ncof)

    # If it exists, extract that information and ask for user confirmation to use it
    if len(act_site) > 0 and use_uniprot == 'YES':
        up_maxid = record.annotations["accessions"][0]
        for i in range(len(act_site)):
            catalytic_res[record.seq[act_site[i]]] = act_site[i]
        copyfile("query.fasta", "refsequp.fasta")
        with open("refsequp.fasta", "a") as f:
            f.write("\n>" + str(record.annotations["accessions"][0])+"\n" + str(record.seq))

    if len(act_site) == 0 or use_uniprot == "NO":
        print("Your model protein does not have an active site identified in UniProt. "
              "Proceeding to check in the M-CSA database...")
        time.sleep(3)
        # Create dictionary from database with uniprot id and EC number only of those
        # sharing the same ECnumber with the model#
        upidEC = {}
        curated_data = csv.reader(open("./curated_data.csv", "r"))
        for row in curated_data:
            upidEC[row[3]] = str(row[1])
        searchUP_dict = {}
        while len(searchUP_dict) == 0:
            for key, value in upidEC.items():
                if proteinECnumber in key:
                    searchUP_dict[key] = value
            else:
                proteinECnumberlist = proteinECnumber.split(".")
                maxi = len(proteinECnumberlist)
                proteinECnumber = str(proteinECnumberlist[0]) + "." + str(proteinECnumberlist[1])
                if proteinECnumber in key:
                    searchUP_dict[key] = value

        # Clean uniprot ids if they have more than 6 characters
        uniprotids = list(searchUP_dict.values())
        uniprotids_clean = []
        for i in uniprotids:
            if len(i) != 6:
                uniprotids_clean.append(i[0:6])
            else:
                uniprotids_clean.append(i)
        # Retrieve sequences of those uniprot identifiers from the website
        upECrefseq = []
        i = 0
        for x in uniprotids_clean:
            handle2 = urllib.request.urlopen("https://www.uniprot.org/uniprot/" + uniprotids_clean[i] + ".xml")
            record2 = SeqIO.read(handle2, "uniprot-xml")
            upECrefseq.append(record2.seq)
            i += 1
        # Use Pairwise alignment to identify the best hit from the M-CSA list
        # Create the fasta file for the identified hit sin the M-CSA
        i = 0
        for seqs in uniprotids_clean:
            with open(uniprotids_clean[i]+".fasta", "w") as f:
                f.write("\n>"+uniprotids_clean[i]+"\n" + str(upECrefseq[i]))
            i += 1

        # Pairwise alignment execution. The options used are the typical from the BLOSUM62 substitution matrix

        aligner = Align.PairwiseAligner()
        alphabet = "PROTEIN"
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.1
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        # Save into a list the scores of each alignment

        alignment_scores = []
        i = 0
        for seqs in uniprotids_clean:
            seq1 = SeqIO.read("query.fasta", "fasta")
            seq2 = SeqIO.read(uniprotids_clean[i]+".fasta", "fasta")
            alignments = aligner.align(seq1.seq, seq2.seq)
            alignment_scores.append(alignments.score)
            i += 1

        # Idenfity the UniProt ID with the maximal identity to the query sequence

        try:
            up_maxid = uniprotids_clean[alignment_scores.index(max(alignment_scores))]
            print("The best hit for your protein is uniprotID: " + str(up_maxid) + ". \n")
            time.sleep(3)
        except BaseException:
            print("Could not find a protein with enough homology with your query. "
                  "You can try again with a different model!")
            time.sleep(4)
            quit()
        time.sleep(2)

        # From the UniProtID with max identity, retrieve the catalytic residues annotated in the M-CSA csv file
        curated_data = csv.reader(open("./curated_data.csv", "r"))
        for row in curated_data:
            if str(row[1])[0:6] == str(up_maxid) and row[4] == "residue":
                if row[5] not in catalytic_res.keys():
                    catalytic_res[row[5]] = [row[7]]
                elif row[5] in catalytic_res.keys() and row[7] not in catalytic_res[row[5]]:
                    catalytic_res[row[5]].append(row[7])

        # From the UniProtID with max identity, retrieve the cofactors (if any)
        curated_data = csv.reader(open("./curated_data.csv", "r"))
        for row in curated_data:
            if str(row[1]) == str(up_maxid) and row[4] == "cofactor":
                cofactor.append(row[8] + "(" + row[5] + ")")
        cofactor = set(cofactor)

    # Alignment of the query sequence with the best hit from UniProt to identify the catalytic residues in the query.
    # Create a fata file containing the query sequence and the sequence of the best UniProtID
        shutil.copyfile("query.fasta", "refsequp.fasta")

        with open("refsequp.fasta", "a") as f:
            f.write("\n>"+up_maxid+"\n" + str(upECrefseq[alignment_scores.index(max(alignment_scores))]))

    # Use ClustalW for the alignment
    sequp_align = ClustalwCommandline(clustalw_exe, infile="refsequp.fasta", score="percent")
    stdout, stderr = sequp_align()
    with open("alignment_output.txt", "w+") as clustalscore:
        clustalscore.write(stdout)

    # From the ClustalW results, check if the identity is sufficient to continue with the active site identification
    with open("alignment_output.txt") as c:
        for line in c:
            if "Sequences (1:2)" in line:
                score = int(''.join(filter(str.isdigit, line)))
    if len(str(score)) == 3:
        score = score - 120
    elif len(str(score)) == 4:
        score = score - 1200
    elif len(str(score)) == 5:
        score = 100

    print("\t Percentage of identity = " + str(score) + "\n")
    time.sleep(3)
    if 40 > score > 15:
        identity = 1
    elif score <= 15:
        identity = 2
    else:
        identity = 0

    time.sleep(3)
    if identity == 0:
        highhomology = "YES"
    elif identity == 1:
        highhomology = "YES"
    else:
        highhomology = "NO"
        print("Your query sequence alignment resulted in less than 15% identity. "
              "Prediction of the active site would not be accurate using it.\n")
        time.sleep(3)

    # From alignment, and the catalytic residues identified in the M-CSA hit or the UniProt,
    # return positions in your protein#
    res = open("Active_site.txt", "a")
    sequp_alignment = AlignIO.read("refsequp.aln", "clustal")

    if highhomology == "YES":
        # Check that the list is not empty
        aaposition_cat = []
        aaposition_cat1 = []
        if len(list(catalytic_res.keys())[0]) >= 1:
            aaname_cat = list(catalytic_res.keys())
            aaposition_cat = list(catalytic_res.values())
            for aaposition in aaposition_cat:
                if isinstance(aaposition, list):
                    for aa in aaposition:
                        aaposition_cat1.append(aa)
                else:
                    aaposition_cat1.append(aaposition)
        
            # Convert the identified residues to a dictionary with the information in a more readable format
            sequp_alignment = AlignIO.read("refsequp.aln", "clustal")
            threetoone = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q",
                          "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
                          "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
            catalytic_res1lett = {}
            i = 0
            conv3to1 = 0
            for keys in catalytic_res.keys():
                if keys in threetoone.keys():
                    conv3to1 = 1
                    catalytic_res1lett[threetoone[str(aaname_cat[i])]] = aaposition_cat[i]
                    i += 1
            if conv3to1 == 1:
                aaname_cat = list(catalytic_res1lett.keys())
            else:
                for positions in catalytic_res.keys():
                    catalytic_res1lett[positions] = [catalytic_res[positions]]  
        else:
            aaname_cat = list(catalytic_res.keys())
            aaposition_cat = list(catalytic_res.values())
            aaposition_cat1 = list(catalytic_res.values())
            for positions in catalytic_res.keys():
                catalytic_res1lett[positions] = [catalytic_res[positions]]
        # From the known catalytic positions, identify the corresponding positions in the alignment
        print("Identifying the predicted catalytic residues in your protein...\n")
        time.sleep(3)
        a = 0
        positions_tocheck = []
        for it in aaposition_cat1:
            positions_tocheck.append(int(it)-1)
            a += 1
        i = 0
        align_catpos = []
        for i in range(a):
            pos_cat = positions_tocheck[i]
            gap_count, gap_count_1, gap_count_2, gap_count_3 = (0, 0, 0, 0)
            for aa in sequp_alignment[1][:positions_tocheck[i]]:
                if aa == "-":
                    gap_count += 1
            if gap_count != 0:
                for aa1 in sequp_alignment[1][pos_cat:pos_cat + gap_count]:   
                    if aa == "-":
                        gap_count_1 += 1
                if gap_count_1 != 0:
                    for aa2 in sequp_alignment[1][pos_cat:pos_cat + gap_count + gap_count_1]:   
                        if aa == "-":
                            gap_count_2 += 1
                    if gap_count_2 != 0:
                        for aa3 in sequp_alignment[1][pos_cat:pos_cat + gap_count + gap_count_1 + gap_count_2]:   
                            if aa == "-":
                                gap_count_3 += 1
            align_catpos.append(pos_cat + gap_count + gap_count_1 + gap_count_2 + gap_count_3)
            i += 1
        # Extract in list from the dictionary of 1lettcode
        lettkeys = list(catalytic_res1lett.keys())

        # From the positions in the alignment, identify the real position in the query sequence
        query_catalytic_res = {}
        for lett in lettkeys:
            query_catalytic_res[lett] = []

        query = SeqIO.read("query.fasta", "fasta")
        for pos in align_catpos:
            pos_cat = pos
            pos_cat_real = pos_cat
            for aas in sequp_alignment[0][0:pos_cat+1]:
                if aas == "-":
                    pos_cat_real -= 1
            try:
                if isinstance(query_catalytic_res[query[pos_cat_real]], list) is True:
                    query_catalytic_res[query[pos_cat_real]].append(pos_cat_real+1)
                elif isinstance(query_catalytic_res[query[pos_cat_real]], str) is True:
                    query_catalytic_res[query[pos_cat_real]] = [query_catalytic_res[query[pos_cat_real]]]
                    query_catalytic_res[query[pos_cat_real]].append(pos_cat_real+1)
            except KeyError:
                try:
                    if isinstance(query_catalytic_res[query[pos_cat_real+1]], list) is True:
                        query_catalytic_res[query[pos_cat_real+1]].append(pos_cat_real+2)
                    elif isinstance(query_catalytic_res[query[pos_cat_real+1]], str) is True:
                        query_catalytic_res[query[pos_cat_real+1]] = [query_catalytic_res[query[pos_cat_real+1]]]
                        query_catalytic_res[query[pos_cat_real+1]].append(pos_cat_real+2)
                except KeyError:
                    try:
                        if isinstance(query_catalytic_res[query[pos_cat_real-1]], list) is True:
                            query_catalytic_res[query[pos_cat_real-1]].append(pos_cat_real)
                        elif isinstance(query_catalytic_res[query[pos_cat_real-1]], str) is True:
                            query_catalytic_res[query[pos_cat_real-1]] = [query_catalytic_res[query[pos_cat_real-1]]]
                            query_catalytic_res[query[pos_cat_real-1]].append(pos_cat_real)
                    except KeyError:
                        query_catalytic_res[query[pos_cat_real]] = str(pos_cat_real+1) + "*"    
    
        query_catalytic_resC = query_catalytic_res.copy()
        for key in query_catalytic_resC.keys():
            if len(query_catalytic_resC[key]) == 0:
                query_catalytic_res[key] = "Present in model and absent in the query."
        time.sleep(3)

        # Print the results both in the console and in a file
        print("The predicted active site is formed by:")
        print("The predicted active site is formed by:", file=res)
        for key, value in query_catalytic_res.items():
            print(str(key) + "--> " + str(value).replace("[", "").replace("]", ""))
            print(str(key) + "--> " + str(value).replace("[", "").replace("]", ""), file=res)
        print("\n* : Residue predicted in the query but not present as part of the active site of the model.\n")
        print("\n* : Residue predicted in the query but not present as part of the active site of the model.\n",
              file=res)
        time.sleep(3)
    else:
        if len(bind_site) > 0:
            print("Due to low homology, calculation of the active site failed. "
                  "Maybe you could check the model pdb file information in  "
                  "https://www.uniprot.org/uniprot/"+str(uniprotID) + " for more information.")
            print("Due to low homology, calculation of the active site failed. "
                  "Maybe you could check the model pdb file information in  "
                  "https://www.uniprot.org/uniprot/"+str(uniprotID) + " for more information.", file=res)
            time.sleep(3)
        else:
            print("Due to low homology, calculation of the active site failed. "
                  "Maybe you could check the model pdb file information in  "
                  "https://www.uniprot.org/uniprot/"+str(uniprotID) + " for more information.")
            print("Due to low homology, calculation of the active site failed. "
                  "Maybe you could check the model pdb file information in  "
                  "https://www.uniprot.org/uniprot/"+str(uniprotID) + " for more information.", file=res)

    # At last, from the EC number write some general information of the protein. If they are supposed to have a cofactor
    # but could not be identified, print also that!
    if len(cofactor) > 0:
        print("\nIdentified cofactor(s): ")
        print("\nIdentified cofactor(s): ", file=res)
        for x in cofactor:
            print(x)
            print(x, file=res)
    elif len(bind_site) > 0:
        print("Also, in Uniprot some binding sites are indicated for the pdb model in positions: ")
        print(*bind_site, sep=",")
        print("Also, in Uniprot some binding sites are indicated for the pdb model in positions: ", file=res)
        print(*bind_site, sep=",", file=res)
    else:
        if "1." in proteinECnumber[0:2]:
            print("Your protein is an oxidoreductase, so it should have a cofactor. But it could not be identified. "
                  "Please, check uniprot entry " + str(up_maxid) + " -the best hit identified from the M-CSA database- "
                                                                   "for more information.")
            print("Your protein is an oxidoreductase, so it should have a cofactor. But it could not be identified. "
                  "Please, check uniprot entry " + str(up_maxid) + " -the best hit identified from the M-CSA database- "
                                                                   "for more information.", file=res)
        elif "2." in proteinECnumber[0:2]:
            print("Your protein is a transferase, commonly, cofactor dependant. Please, check uniprot entry "
                  + str(up_maxid) + " -the best hit identified- for more information.")
            print("Your protein is a transferase, commonly, cofactor dependant. Please, check uniprot entry "
                  + str(up_maxid) + " -the best hit identified- for more information.", file=res)
        elif "3." in proteinECnumber[0:2]:
            print("Your protein is an hydrolase.")
            print("Your protein is an hydrolase.", file=res)
        elif "4." in proteinECnumber[0:2]:
            print("Your protein is most similar to characterised lyases.")
            print("Your protein is most similar to characterised lyases.", file=res)
        elif "5." in proteinECnumber[0:2]:
            print("Your protein seems to be an isomerase.")
            print("Your protein seems to be an isomerase.", file=res)
        elif "6." in proteinECnumber[0:2]:
            print("Your protein seems to be an ligases. Normally, this enzymes are ATP dependant but no cofactor could "
                  "be identified. Please, check uniprot entry " + str(up_maxid) +
                  " -the best hit identified- for more information.")
            print("Your protein seems to be an ligases. Normally, this enzymes are ATP dependant but no cofactor could "
                  "be identified. Please, check uniprot entry " + str(up_maxid) +
                  " -the best hit identified- for more information.", file=res)

    res.close()

    time.sleep(3)

    # Clean the directory of intermediate files which are not necessary
    os.remove("refsequp.dnd")
    tokeep = ["query.fasta", str(up_maxid)+".fasta"]
    for file in glob.glob("*.fasta"):
        if file not in tokeep:
            os.remove(file)

    print("\nFinished running the ActiveSiteID!\n. Your results can be found in " + str(os.getcwd()) +
          ", in a file named Active_site.txt, together with the alignment and the fasta files for your sequences. "
          "You can run now the Surface&Clusters module! ")

    os.chdir(initial_location)
   
    return uniprotID

if __name__ == "__main__":
    main(folder_name, folder_answer, prot_sequence, model_prot, use_uniprot)
