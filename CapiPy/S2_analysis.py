#!/usr/bin/env python3

# Size&Clusters v.0.1
# Copyleft David Roura, 2020
"""Surface analysis and residue clustering of a query protein based on it's three dimensional structure
"""

import os
import shutil
import time
import csv
import numpy as np
import re
import random

from os.path import dirname
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *
from modeller import *
from modeller.automodel import *
from random import randint
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)


def main(folder_answer, folder_name, pdbfile, active_site):
    initial_location = os.getcwd()
    # Locate the workspace
    working_directory = folder_name
    if folder_answer == 'NO':
        if os.path.exists(working_directory) is True:
            os.chdir(working_directory)
            if os.path.exists("Size&Clusters") is False:
                os.mkdir("Size&Clusters")
            os.chdir("./Size&Clusters")
        try:
            shutil.copy("../Blast&Modeller/" + str(pdbfile), './query.pdb')
        except BaseException:
            print('Cannot find the pdb file. Check again the name!')
            quit()
    elif folder_answer == 'YES':
        os.mkdir(working_directory)
        os.chdir(working_directory)
        os.mkdir("Size&Clusters")
        os.chdir("./Size&Clusters")
        try:
            shutil.copy(pdbfile, './query.pdb')
        except BaseException:
            print('Cannot find the pdb file. Check again the name!')
            quit()

    # Open the desired pdb file
    parser = PDBParser()
    fullid = "Q00F"
    fullfile = "query.pdb"
    full_structure = parser.get_structure(fullid, fullfile)
    mdel = full_structure[0]
    ppb = PPBuilder()

    # Calculate the volume of the minimal bounding box
    ca_coord = []
    number_of_chains = 0
    for chain in mdel.get_list():
        number_of_chains += 1
        for residue in chain.get_list():
            if residue.has_id("CA"):
                ca = residue["CA"]
                ca_coord.append(ca.get_coord())

    print("\nMeasuring the minimum bounding box...\n")

    ca_coord = np.array(ca_coord)
    min_x = min(ca_coord[:, 0])
    min_y = min(ca_coord[:, 1])
    min_z = min(ca_coord[:, 2])
    max_x = max(ca_coord[:, 0])
    max_y = max(ca_coord[:, 1])
    max_z = max(ca_coord[:, 2])
    box_dimensions = ([max_x - min_x, max_y - min_y, max_z - min_z])
    time.sleep(2)
    volume = box_dimensions[0] * box_dimensions[1] * box_dimensions[2]
    print("Your protein, formed by " + str(number_of_chains) + " chains, has an aproximate dimensions of: "
          + str(box_dimensions[0]) + " A, " + str(box_dimensions[1]) + " A, " + str(box_dimensions[2]) +
          " A.\nGiving an aproximate volume of: " + str(volume) + " A\u00b3\n")
    time.sleep(3)
    # Sort the amino acids of the structure depending on their exposure calculated by the half-sphere exposure method
    # - 10.1093/bioinformatics/btv665
    print("Calculating the exposure of each residue in the structure...")
    RADIUS = 16.0
    RADIUS = 16.0
    hse_cb = HSExposureCB(mdel, radius=RADIUS)
    buried = []
    semiburied = []
    exposed = []
    print("Sorting the residues depending on their exposure...")
    for r in mdel.get_residues():
        if is_aa(r) and r.xtra["EXP_HSE_B_U"] >= 35 and r.xtra["EXP_HSE_B_D"] >= 35:
            buried.append(r)
        elif is_aa(r) and 35 > r.xtra["EXP_HSE_B_U"] >= 25 and 35 > r.xtra["EXP_HSE_B_D"] >= 25:
            semiburied.append(r)
        elif is_aa(r):
            exposed.append(r)

    # Calculate the surface in A according to the predicted surface of each amino acid classified as exposed
    # - 10.1371/journal.pone.0080635
    aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
          "SER", "THR", "TRP", "TYR", "VAL"]
    aa_s2 = {}
    aa_n = {}
    for aaname in aa:
        aa_s2.update({aaname: int()})
        aa_n.update({aaname: int()})

    aa_ASA = {"ALA": 121, "ARG": 265, "ASN": 187, "ASP": 187, "CYS": 148, "GLU": 214, "GLN": 214, "GLY": 97,
              "HIS": 216, "ILE": 195, "LEU": 191, "LYS": 230, "MET": 203, "PHE": 228, "PRO": 154, "SER": 143,
              "THR": 163, "TRP": 265, "TYR": 255, "VAL": 165}
    totals2 = 0
    for r in exposed:
        totals2 += aa_ASA[r.get_resname()]
        i = 0
        for i in range(20):
            if r.get_resname() == aa[i]:
                aa_s2[aa[i]] += aa_ASA[aa[i]]
                aa_n[aa[i]] += 1
            else:
                i += 1

    aa_s2_perc = {}
    for aaname in aa:
        aa_s2_perc.update({aaname: str("{:.2f}".format(((aa_s2[aaname] / totals2) * 100)) + "%")})

    time.sleep(4)

    print("\nAs for the surface, the calculated exposed surface is " + str(totals2) + " A\u00b2\n")

    print("The contribution of each aminoacid is:")
    print("Aminoacid\t Surface(%)\t Total number")
    for aa in aa_s2_perc:
        print(aa + "\t\t" + aa_s2_perc[aa] + "\t\t" + str(aa_n[aa]))

    print("Total number of surface residues: " + str(len(exposed)) + " out of " + str(len(exposed) + len(semiburied) + len(buried)))

    time.sleep(5)
    # Check if the files already exists, not to have a repetition
    if os.path.exists("General_info.txt"):
        check = open("General_info.txt").read()
    else:
        check = ""
    if "As for" not in check:
        with open("General_info.txt", "a") as f1:
            print("Your protein aproximate dimensions in Angstroms are: " + str(box_dimensions[0]) + "A, " +
                  str(box_dimensions[1]) + "A, " + str(
                box_dimensions[2]) + "A.\nWhich makes an aproximate volume of: " +
                  str(volume) + " A\u00b3\\n", file=f1)
            print("\nAs for its surface, it measures " + str(totals2) + " A\u00b2\\n", file=f1)
            print("From which each amino acid contributes to:", file=f1)
            print("Aminoacid\t Surface(%)\t Total number", file=f1)
            for aa in aa_s2_perc:
                print(aa + "\t\t" + aa_s2_perc[aa] + "\t\t" + str(aa_n[aa]), file=f1)
            print("Total number of surface residues: " + str(len(exposed)), file=f1)

    time.sleep(3)

    # Clustering part: First, sort amino acids into groups of the most important types for protein immobilisation;
    # Second, check if those amino acids are at less than 10A. If so, identify them as a cluster.
    aatocheck = ["LYS", "GLU", "ASP", "HIS", "CYS", "TYR", "ARG"]  # commonly used for immobilisation#
    hydrophobic = ["ILE", "LEU", "PHE", "VAL", "TRP", "TYR"]  # arunachalam 2008#
    charged = []
    hydroph = []

    for r in exposed:
        name = r.get_resname()
        if name in aatocheck:
            chainid = r.get_parent().get_id()
            resseq = r.get_id()[1]
            charged.append(name + "_" + str(resseq) + "_" + chainid)
        elif name in hydrophobic:
            chainid = r.get_parent().get_id()
            resseq = r.get_id()[1]
            hydroph.append(name + "_" + str(resseq) + "_" + chainid)

    poslist = []
    neglist = []
    hislist = []
    lyslist = []
    cyslist = []

    for item in charged:
        if "LYS" in item:
            poslist.append(item)
            lyslist.append(item)
        elif "ASP" in item or "GLU" in item:
            neglist.append(item)
        elif "HIS" in item:
            poslist.append(item)
            hislist.append(item)
        elif "CYS" in item:
            cyslist.append(item)
        elif "ARG" in item:
            poslist.append(item)

    # ______________________________CLUSTERING_________________________________________________________________________

    def clustering(list1, dic):
        # This method checks if the residues in a certain list have their CA at less than 10A.
        # If so, this residues will be included in the dictionary as a new cluster
        full_structure = parser.get_structure(fullid, fullfile)
        i = 0
        while i in range(len(list1)):
            x1 = list1[i]
            for x2 in list1:
                if x1 != x2:
                    dis = mdel[x1.split("_")[2]][int(x1.split("_")[1])]["CA"] - \
                          mdel[x2.split("_")[2]][int(x2.split("_")[1])]["CA"]
                    if dis < 10:
                        try:
                            dic[x1].append(x2)
                        except KeyError:
                            dic[x1] = [x2]
            i += 1

    def clean(dict1):
        # Clean up the dictionaries created with clustering. This avoids repetition of the same cluster as well as
        # deletes any cluster formed by less than 3 different aminoacids
        dictcopy = dict1.copy()
        keys = list(dictcopy.keys())
        for key in keys:
            if key in dict1[key]:
                dict1.pop(key)
            elif len(dict1[key]) < 2:
                dict1.pop(key)
            elif key in str(dict1.values()):
                dict1.pop(key)
            elif "" in dict1[key]:
                while "" in dict1[key]:
                    dict1[key].remove("")

    def writeclusterstocsv(cluster, file):
        # Self explanatory. Writes the clusters identified into a csv file.
        with open(file, "a", newline='') as f:
            writer = csv.writer(f)
            if len(cluster) > 1:
                for i in cluster:
                    towrite = [i]
                    for x in cluster[i]:
                        towrite.append(x)
                    writer.writerow(towrite)
                writer.writerow("\n")

    def show_clusters(dict):
        # Writes the information of teh clusters into a txt file. More "user friendly" explanation of the results
        f = open("Clusters.txt", "a")
        if len(dict) > 1:
            f.write(list(dict.keys())[0] + " residues in a cluster (<10A):\n")
            i = 1
            for keys in list(dict.keys())[1:]:
                toprint = "Cluster " + str(list(dict.keys())[0][0:3].lower()) + str(i) + ":" + str(keys)
                length = len(dict[keys])
                x = 0
                while x < length:
                    toprint += ", " + str(dict[keys][x])
                    x += 1
                toprint += "."
                print(toprint, file=f)
                i += 1
        f.write("\n")

    def uploadcluster():
        # In case the clusters.csv (which is created in the first execution -  is available, file is opened and clusters
        # written back into a Python dictionary to continue execution
        csv_clusters = csv.reader(open("clusters.csv", "r"))
        i = 0
        keys = [[], [], [], [], [], []]
        values = [[], [], [], [], [], []]
        with open("clusters.csv", "r") as file:
            csv_clusters = csv.reader(file)
            for row in csv_clusters:
                if row[0] == "\n":
                    i += 1
                else:
                    keys[i].append(row[0])
                    values[i].append(row[1:])
        x = 0
        for x in range(len(keys)):
            try:
                if keys[x][0] == "Positive":
                    for y in range(len(keys[x])):
                        cluster_pos[keys[0][y]] = values[x][y]
                elif keys[x][0] == "Negative":
                    for y in range(len(keys[x])):
                        cluster_neg[keys[x][y]] = values[x][y]
                elif keys[x][0] == "Histidine":
                    for y in range(len(keys[x])):
                        cluster_his[keys[x][y]] = values[x][y]
                elif keys[x][0] == "Lysine":
                    for y in range(len(keys[x])):
                        cluster_lys[keys[x][y]] = values[x][y]
                elif keys[x][0] == "Cysteine":
                    for y in range(len(keys[x])):
                        cluster_cys[keys[x][y]] = values[x][y]
                elif keys[x][0] == "Hydrophobic":
                    for y in range(len(keys[x])):
                        cluster_hidroph[keys[x][y]] = values[x][y]
            except IndexError:
                x += 1
            x += 1

    def clusters_pymol(dict):
        # Writes the clusters into a pymol command which can be direclty copied into the command line
        if len(dict) > 1:
            f = open("PyMol_clusters.pml", "a")
            colours = coloursdic[list(dict.keys())[0]]
            i = 0
            for keys in dict.keys():
                if keys.isalpha() is not True:
                    toprint = ("select " + str(list(dict.keys())[0][0:3].lower()) + str(i) + ", (resi " + str(
                        keys.split("_")[1])) + " & chain " + str(dict[keys][1].split("_")[2])
                    length = len(dict[keys])
                    x = 0
                    while x < length:
                        toprint += ", resi " + str(dict[keys][x].split("_")[1]) + " & chain " + str(
                            dict[keys][x].split("_")[2])
                        x += 1
                    toprint += ")"
                    print(toprint, file=f)
                i += 1
            print("color " + colours + ", " + str(list(dict.keys())[0][0:3].lower()) + "*\n", file=f)

    # Start with the calculation of the clusters, if user agrees
    cluster_pos = {"Positive": ["Lys", "Arg", "His"]}
    cluster_neg = {"Negative": ["Asp", "Gly"]}
    cluster_his = {"Histidine": ["only", "histidines"]}
    cluster_lys = {"Lysine": ["only", "lysines"]}
    cluster_cys = {"Cysteine": ["only", "cysteines"]}
    cluster_hidroph = {"Hydrophobic": ["ILE", "LEU", "PHE", "VAL", "TRP", "TYR"]}
    coloursdic = {"Positive": "Blue", "Negative": "Red", "Histidine": "Cyan", "Lysine": "LightBlue",
                  "Cysteine": "Yellow", "Hydrophobic": "White"}

    if os.path.exists("clusters.csv") is True:
        print("Uploading previously calculated clusters form clusters.csv.\n")
        uploadcluster()
    elif os.path.exists("clusters.csv") is False:
        print("Measuring distances and identifying clusters...")
        clustering(poslist, cluster_pos)
        clustering(neglist, cluster_neg)
        clustering(hislist, cluster_his)
        clustering(lyslist, cluster_lys)
        clustering(cyslist, cluster_cys)
        clustering(hydroph, cluster_hidroph)
        clean(cluster_pos)
        clean(cluster_neg)
        clean(cluster_his)
        clean(cluster_lys)
        clean(cluster_cys)
        clean(cluster_hidroph)
        writeclusterstocsv(cluster_pos, "clusters.csv")
        writeclusterstocsv(cluster_neg, "clusters.csv")
        writeclusterstocsv(cluster_his, "clusters.csv")
        writeclusterstocsv(cluster_lys, "clusters.csv")
        writeclusterstocsv(cluster_cys, "clusters.csv")
        writeclusterstocsv(cluster_hidroph, "clusters.csv")
        show_clusters(cluster_pos)
        show_clusters(cluster_neg)
        show_clusters(cluster_his)
        show_clusters(cluster_lys)
        show_clusters(cluster_cys)
        show_clusters(cluster_hidroph)
        print("Clusters can be found in Clusters.txt")
    time.sleep(2)
    with open("PyMol_clusters.pml", "w") as f:
        print("load model.pdb", file=f)
        chainid = []
        for chains in full_structure.get_chains():
            if chains.get_id() not in chainid:
                chainid.append(chains.get_id())
        if len(chains) > 1:
            for chain in chainid:
                print("set_color shade" + str(chain) + "= [1.0, " + str(random.random()) + ", " + str(
                    random.random()) + "]", file=f)
                print("colour shade" + str(chain) + ", chain " + str(chain), file=f)
    clusters_pymol(cluster_pos)
    clusters_pymol(cluster_neg)
    clusters_pymol(cluster_his)
    clusters_pymol(cluster_lys)
    clusters_pymol(cluster_cys)
    clusters_pymol(cluster_hidroph)
    time.sleep(2)
    print("\nClusters can be visualised in PyMOL opening PyMol_clusters.pml\n")
    print("Default colours are:")
    for i, j in coloursdic.items():
        print(i + ": " + j)
    time.sleep(3)

    # ____________________ MULTIMERIC INTERFACE + CLUSTERS __________________________________________________________

    def distointer(cluster, list1):
        # Check if the residues in the cluser can interfere with the protein interface
        itemslist = []
        for kcluster in list(cluster.keys())[1:]:
            itemslist.append(kcluster)
            for values in cluster[kcluster]:
                itemslist.append(values)
        for pos1 in itemslist:
            for pos2 in interface_res["A"]:
                dis = mdel[str(pos1.split("_")[2])][int(pos1.split("_")[1])]["CA"] - mdel["A"][pos2]["CA"]
                if dis != 0 and dis <= 10 and (pos1 not in list1):
                    list1.append(pos1)
                if IndexError:
                    continue

    def writecluster(original, d1, list1, file):
        # Similar to show_cluster method, it writes the information into a .txt file. In this case, it will add a
        # -Warning- if the cluster is close to the specified residues.'''
        it = 1
        f = open(file, "a")
        with open(file) as readfile:
            if "Chain" not in readfile.read():
                f.write("Searching for clusters in close contact to: ")
                if isinstance(original, dict) is True:
                    for aas in original:
                        f.write("\nChain " + str(aas) + ":")
                        for aas2 in original[aas]:
                            f.write(" - " + str(aas2))
                elif isinstance(original, list):
                    for aas in original:
                        f.write(" - " + str(aas))
                readfile.close()
            else:
                readfile.close()
        f.write("\n" + list(d1.keys())[0] + " residues in a cluster (<10A):\n")
        inlist1 = []
        for k in list(d1.keys())[1:]:
            toprint = "Cluster " + str(list(d1.keys())[0][0:3].lower()) + str(it) + ":" + str(k)
            length = len(d1[k])
            inlist1.append(k)
            x = 0
            while x < length:
                toprint += ", " + str(d1[k][x])
                inlist1.append(d1[k][x])
                x += 1
            toprint += "."
            for ex in inlist1:
                if ex in list1 and "\t -In close proximity to specified residue/s - " not in toprint:
                    toprint += "\t -In close proximity to specified residue/s - "
            print(toprint, file=f)
            it += 1
        f.write("\n")
        f.close()

    # Start first, identifying how many chains does the model have. Maximum of 10 for convenience
    # (and because most proteins fall into this category).
    chainres = {"A": [], "B": [], "C": [], "D": [], "E": [], "F": [], "G": [], "H": [], "I": [], "J": []}
    chainids = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]

    # number of chains in model
    number_of_chains = 0
    for chain in mdel.get_chains():
        number_of_chains += 1

    # Create dictionary with all CA location of a chain
    for x in range(number_of_chains):
        for residues in mdel[chainids[x]]:
            if is_aa(residues) is True:
                chainres[chainids[x]].append(residues["CA"])

    # Delete unnecessary keys in the chainres dictionary
    for a in range(len(chainres.keys())):
        if len(chainres[chainids[a]]) == 0:
            del chainres[chainids[a]]

    # calculation of the interface residues. This will assume all chains will have the same interactions as chain A.
    # This simplification is necessary to get results in a reasonable time.
    interface_res = {"A": []}

    if number_of_chains > 1:
        print("Calculating if any of the clusters might affect the quaternary structure of your protein. \n "
              "Note: This assumes chain A is in contact with all other subunits and all contacts are identical "
              "between subunits.\n")

    calculationdone = os.path.exists("residues_interface.txt")

    print("Calculating contacts between chain A and it's neighbours...")
    ch_1 = 1
    # Calculate the distance of all atoms in chain A to all atoms in other chains. If this distance is <10A,
    # the position will be added as part of the interface.
    while ch_1 < number_of_chains:
        for ca1 in chainres["A"]:
            for ca2 in chainres[chainids[ch_1]]:
                if ca1.get_parent().get_id()[1] not in interface_res["A"]:
                    dist = ca1 - ca2
                    if 0 < dist < 10:
                        interface_res["A"].append(ca1.get_parent().get_id()[1])
        ch_1 += 1

        # Assuming all contacts in the different chains are the same as in chain A,
        # copy the residue position to each chain.
        count_start = 1
        for count_start in range(number_of_chains):
            interface_res[chainids[count_start]] = interface_res["A"]
            count_start += 1

        # Write information to file so it doesn't have to be calculated again.
        with open("residues_interface.txt", "a") as fileint:
            print(interface_res, file=fileint)
        print("The position of the residues identified as part of the interface can be found in residues_interface.txt")

    # Calculate now if any of the identified clusters is at less than 10A to the any of the interface residues.
    intpos = []
    intneg = []
    inthis = []
    intlys = []
    intcys = []
    inthidroph = []

    distointer(cluster_pos, intpos)
    writecluster(interface_res, cluster_pos, intpos, "Interface_contacts.txt")
    distointer(cluster_neg, intneg)
    writecluster(interface_res, cluster_neg, intneg, "Interface_contacts.txt")
    distointer(cluster_his, inthis)
    writecluster(interface_res, cluster_his, inthis, "Interface_contacts.txt")
    distointer(cluster_lys, intlys)
    writecluster(interface_res, cluster_lys, intlys, "Interface_contacts.txt")
    distointer(cluster_cys, intcys)
    writecluster(interface_res, cluster_cys, intcys, "Interface_contacts.txt")
    distointer(cluster_hidroph, inthidroph)
    writecluster(interface_res, cluster_hidroph, inthidroph, "Interface_contacts.txt")
    print("Distance to the multimeric interface successful. Your results can be found in Cluster_interface.txt")

    # __________________________ACTIVE SITE + CLUSTERS _________________________________________________________________
    #   
    def dist_to_as(dict1, list1):
        # Similarly to the other distance methods, calculates if any residues of the cluster is at less than
        # 10A of any of the residues defined as the active site'''
        cluster_aa = []
        for k in list(dict1.keys())[1:]:
            cluster_aa.append(k)
            for i in dict1[k]:
                cluster_aa.append(i)
        for as_site in aspositions:
            for x in cluster_aa:
                for i in range(number_of_chains):
                    dis = mdel[chainids[i]][int(as_site)]["CA"] - mdel[x.split("_")[2]][int(x.split("_")[1])]["CA"]
                    if dis < 5:
                        if x not in list1:
                            list1.append(x)

    if active_site == "YES":
        # If this information exists in the active_site directory, give the option to use directly that information.
        # If not, ask for user input on where the active site is located
        parent_dir = dirname(os.getcwd().replace("\\", "/"))
        str_ac = parent_dir + "/ActiveSite/Active_site.txt"
        aspositions_0 = []
        if os.path.exists(str_ac) is True:
            with open(parent_dir + "/ActiveSite/Active_site.txt", "r") as f:
                for line in f.readlines():
                    if "-->" in line:
                        toadd = re.findall(r'\d+', line)
                        aspositions_0.append(toadd)
        else:
            print("Cannot find the active site identified with the ActiveSiteID module, "
                  "please run it and check the results!")

    else:
        print("Using user defined positions for the active site...")
        aspositions_0 = active_site
        
    aspositions = []
    for aa in aspositions_0:
        if isinstance(aa, list) is True:
            for aa1 in aa:
                aspositions.append(aa1)
        else:
            aspositions.append(aa)

        # Calculate if any residue in the cluster is close to the active site residues
    asposd = []
    asnegd = []
    ashisd = []
    aslysd = []
    ascysd = []
    ashidrophd = []

    number_of_chains = 0
    for chain in mdel.get_chains():
        number_of_chains += 1
        
    dist_to_as(cluster_pos, asposd)
    dist_to_as(cluster_neg, asnegd)
    dist_to_as(cluster_his, ashisd)
    dist_to_as(cluster_lys, aslysd)
    dist_to_as(cluster_cys, ascysd)
    dist_to_as(cluster_hidroph, ashidrophd)
    writecluster(aspositions, cluster_pos, asposd, "Cluster_ActiveSite.txt")
    writecluster(aspositions, cluster_neg, asnegd, "Cluster_ActiveSite.txt")
    writecluster(aspositions, cluster_his, ashisd, "Cluster_ActiveSite.txt")
    writecluster(aspositions, cluster_lys, aslysd, "Cluster_ActiveSite.txt")
    writecluster(aspositions, cluster_cys, ascysd, "Cluster_ActiveSite.txt")
    writecluster(aspositions, cluster_hidroph, ashidrophd, "Cluster_ActiveSite.txt")
    print("Information of the cluster distance to the active site can be found in Cluster_ActiveSite.txt")

    # Clean the directory of intermediate files which are not necessary
    time.sleep(2)

    if os.path.exists("query_multimeric.pdb"):
        os.rename("query_multimeric.pdb", "model.pdb")
    elif os.path.exists("query_monomer.pdb"):
        os.rename("query_monomer.pdb", "model.pdb")

    print("\n\tFinished running the Surface&Clusters module!\nYou can find your result files in " + str(
        os.getcwd()) + "in the folder Size&Clusters. For more information, check the tutorial.\n\t "
                       "OH, and make sure to check out PyMol_clusters.pml!")

    os.chdir(initial_location)


if __name__ == "__main__":
    main(folder_answer, folder_name, pdbfile, active_site)
