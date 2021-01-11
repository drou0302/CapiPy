#!/usr/bin/env python3

# Clustering v.0.1
# Copyleft David Roura, 2020
"""Distance calculator between identified clusters and a user specified residue(s)

"""

import os
import time
import csv
import glob

from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *
from modeller import *
from modeller.automodel import *
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)


# Locate the workspace#


def main(folder_name, residues, file_name):
    initial_location = os.getcwd()
    # Locate the workspace#
    working_directory = folder_name
    try:
        if os.path.exists(working_directory) is True:
            os.chdir(working_directory)
            if os.path.exists("Size&Clusters") is False:
                print("You have to have run the third module, Surface and Clusters, before running this one!")
            else:
                os.chdir("Size&Clusters")
    except OSError:
        print("Error accessing the folder. Are you sure it exists?")

    # Open the desired pdb file
    parser = PDBParser()
    fullid = "Q00F"
    pdbfiles = []
    fullfile = 'query.pdb'
    full_structure = parser.get_structure(fullid, fullfile)
    mdel = full_structure[0]
    ppb = PPBuilder()

    def uploadcluster():
        # Clusters written back into a Python dictionary to continue execution
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

    def distonrl(cluster, list1):
        # Similarly to the other distance methods, calculates if any residues of the cluster is at less than 10A
        # of any of the residues defined by the user
        itemslist = []
        for k0 in list(cluster.keys())[1:]:
            itemslist.append(k0)
            for values in cluster[k0]:
                itemslist.append(values)
        for i1 in itemslist:
            for k1 in nrl_dic.keys():
                for v0 in nrl_dic[k1]:
                    dis = mdel[str(i1.split("_")[2])][int(i1.split("_")[1])]["CA"] - mdel[str(k1)][int(v0)]["CA"]
                    if dis != 0 and dis < 10 and i1 not in list1:
                        list1.append(i1)
                    elif IndexError:
                        continue

    def writecluster(original, d1, list1, file):
        # Writes the information into a .txt file. In this case, it will add a -Warning- if the cluster is close
        # to the specified residues.
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

    def distance_pml(dict):
        n_pmls = []
        for file in glob.glob("*.pml", recursive=True):
            n_pmls.append(file)
        if len(n_pmls)>0:
            filename = "ClusterDistance" + str(len(n_pmls)-1) + ".pml"
        else:
            filename= "ClusterDistance.pml"
        # Writes the clusters into a pymol command which can be directly copied into the command line
        for key in dict.keys():
            writing = 'select res_group' + str(len(n_pmls)) + ', ('
            j=0
            item_n = 0
            for aa_n in dict[key]:
                length = len(dict[key])
                if item_n == 0:
                    writing += 'resi ' + aa_n + ' & chain ' + key
                    item_n += 1
                else:
                    writing += ', resi ' + aa_n + ' & chain ' + key
                    item_n += 1
            writing += ')\n'
            writing += 'zoom (res_group' + str(len(n_pmls)) + ')'
            with open(filename, "a") as f:
                f.write(writing)
        return filename
    # Extract from the pdb file the chain id and the sequence of each chain saved into a dictionary.

    cluster_pos = {"Positive": ["Lys", "Arg", "His"]}
    cluster_neg = {"Negative": ["Asp", "Gly"]}
    cluster_his = {"Histidine": ["only", "histidines"]}
    cluster_lys = {"Lysine": ["only", "lysines"]}
    cluster_cys = {"Cysteine": ["only", "cysteines"]}
    cluster_hidroph = {"Hydrophobic": ["ILE", "LEU", "PHE", "VAL", "TRP", "TYR"]}
    uploadcluster()
    c_ids = {}
    for chains in mdel:
        c_ids[chains.get_id()] = []
    for chains in c_ids:
        c_ids[chains] = [ppb.build_peptides(mdel[chains])[0][0].get_id()[1],
                         ppb.build_peptides(mdel[chains])[0][-1].get_id()[1]]
    # Ask user for the specific residues to check. Check that the specified residues exist in the pdb file.
    # nrl = new residues list, but I wanted to keep it short.
    nrl = residues
    for res in nrl:
        if res.split("_")[1] not in c_ids.keys():
            print("\nIt seems that your input is part of an unexistent chain.")
            quit()
        elif int(res.split("_")[0]) not in range(c_ids[res.split("_")[1]][0], c_ids[res.split("_")[1]][1]):
            print("\nIt seems that your input is not part of the pdb file.")
            quit()
    nrl_dic = {}
    for posit in nrl:
        if posit.split("_")[1] in nrl_dic.keys():
            nrl_dic[posit.split("_")[1]].append(posit.split("_")[0])
        else:
            nrl_dic[posit.split("_")[1]] = [posit.split("_")[0]]
    perspos = []
    persneg = []
    pershis = []
    perslys = []
    perscys = []
    pershidroph = []

    distonrl(cluster_pos, perspos)
    distonrl(cluster_neg, persneg)
    distonrl(cluster_his, pershis)
    distonrl(cluster_lys, perslys)
    distonrl(cluster_cys, perscys)
    distonrl(cluster_hidroph, pershidroph)

    spfile = file_name

    writecluster(nrl_dic, cluster_pos, perspos, spfile + ".txt")
    writecluster(nrl_dic, cluster_neg, persneg, spfile + ".txt")
    writecluster(nrl_dic, cluster_his, pershis, spfile + ".txt")
    writecluster(nrl_dic, cluster_lys, perslys, spfile + ".txt")
    writecluster(nrl_dic, cluster_cys, perscys, spfile + ".txt")
    writecluster(nrl_dic, cluster_hidroph, pershidroph, spfile + ".txt")

    pml_file = distance_pml(nrl_dic)
    print("\nALL DONE, you can find your results files in " + str(os.getcwd()) + " in the folder Size&Clusters")

    os.chdir(initial_location)
    return pml_file

if __name__ == "__main__":
    main(folder_name, residues, file_name)
