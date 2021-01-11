#!/usr/bin/env python3

# Blast&modeller v.0.1
# Copyleft David Roura, 2020
"""Superimposer to mimic a quaternary assembly from a monomeric structure 

"""

import glob
import os
import shutil
import numpy
import time

from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO, SearchIO, PDB, AlignIO
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *
from modeller import *
from modeller.automodel import *
from sys import platform
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)


def main(folder_name, monomeric_model, multimeric_model, clustalw_exe):
    initial_location = os.getcwd()
# Locate the workspace#
    working_directory = folder_name
    try:
        if os.path.exists(working_directory) is True:
            os.chdir(working_directory)
            if os.path.exists("Blast&Modeller") is False:
                os.mkdir("Blast&Modeller")
            os.chdir("./Blast&Modeller")
    except OSError:
        print("Error accessing the folder. Are you sure it exists?")
        quit()
    protein_sequence = monomeric_model
    try:
        shutil.copy(protein_sequence, "./monomeric_structure.pdb")
    except OSError:
        print("Monomeric model not found. Please try to answer again.")
    answer_pdb_id = multimeric_model
    if len(answer_pdb_id) == 4:
        if answer_pdb_id.isalpha() is True:
            print("It seems that your format is not a PDBid. Try again.")
        elif answer_pdb_id.isnumeric() is True:
            print("Remember that PDBid are formed by letters and numbers.")
        else:
            id_1 = answer_pdb_id
            print(str(id_1) + " is going to be used as the template for your model.")
    else:
        print("That doesn't seem to be an accepted format.")

    pdbl = PDBList()
    pdbl.retrieve_pdb_file(answer_pdb_id, pdir=".", file_format="pdb")
    parser = PDBParser()
    ppb = PPBuilder()
    pdbfile = "pdb" + answer_pdb_id + ".ent"
    file_name = str(answer_pdb_id.casefold())
    structure = parser.get_structure(answer_pdb_id, pdbfile)
    env = environ()

    def superimpose_chains(dict):    
        checks = {}
        # Check if clustalw is available.
        for k_chain in toalign.keys():
            chainmodel = parser.get_structure("CH0" + toalign[k_chain], "chain" + toalign[k_chain] + ".pdb")
            chainquery = parser.get_structure("Q00" + str(k_chain), "query_aligned"+str(k_chain) + ".pdb")
            ppb = PPBuilder()
            # To maximise the accuracy of the superimposer and avoid problems with alignments of low quality instead of
            # using the whole protein to superimpose, two regions at the start and end will be used. To do so,
            # first it's necessary to identify two regions at the start and end where both there are no gaps,
            # checking the alignment for a minimum of 40 positions without "-" in reading frames of 5 aminoacids.
            frac_chainmodelseq = []
            chainmodelseq = ""
            for pp in ppb.build_peptides(chainmodel):
                frac_chainmodelseq.append(pp.get_sequence())
            for fractions in frac_chainmodelseq:
                chainmodelseq += fractions
            for pp in ppb.build_peptides(chainquery):
                chainqueryseq = pp.get_sequence()
            with open("initial" + str(k_chain) + ".fasta", "w") as input_file:
                input_file.write(">model_" + toalign[k_chain] + "\n" + str(chainmodelseq) + "\n>query_" + str(k_chain)
                                 + "\n" + str(chainqueryseq))

            ex_clustal = ClustalwCommandline(clustalw_exe, infile="initial"+str(k_chain) + ".fasta")
            stdout, stderr = ex_clustal()
            alignment = AlignIO.read("initial" + str(k_chain) + ".aln", "clustal")
            length = len(alignment[0]) 
            # Identify the first region of 40 aa without gaps in the alignment
            i = 0
            posm = 0
            posq = 0
            posmf = 0
            posqf = 0
            for a in range(int(length/5)):
                if "-" not in alignment[0][i:i+40].seq and "-" not in alignment[1][i:i+40]:
                    posq = i
                    posm = i
                else:
                    i += 5
            # Translate the alignment position to the real position in each sequence
            for aa in alignment[0].seq[0:i]:
                if aa == "-":
                    posm -= 1
            for aa in alignment[1].seq[0:i]:
                if aa == "-":
                    posq -= 1
            # Identify in the last quarter of the the alignment, a region of 40 aa without gaps
            k = length - (length//4)
            while k < length-40:
                if "-" not in alignment[0][k:k+40].seq and "-" not in alignment[1][k:k+40].seq:
                    posqf = k
                    posmf = k
                    break
                else:
                    k += 5
            if k+40 > length:
                k = length - (length//4)
                while k < length-40:
                    if "-" not in alignment[0][k:k+40].seq and "-" not in alignment[1][k:k+40].seq:
                        posqf = k
                        posmf = k
                        break
                    else:
                        k += 2
            if posqf == 0:
                posqf = length - (length//4)
                posmf = length - (length//4)
            # Translate the alignment position to the real position in each sequence
            gapsmodel = 0
            for aa in alignment[0].seq[0:k]:
                if aa == "-":
                    gapsmodel += 1
                    posmf -= 1
            for aa in alignment[1].seq[0:k]:
                if aa == "-":
                    posqf -= 1
            # Check if the pdb file starts at amino acid 1 or not. If not, fix the range to be applied to the model.
            mdel_res_list = PDB.Selection.unfold_entities(chainmodel, "R")
            first_aa_model = mdel_res_list[0].get_id()[1]
            if first_aa_model != 1:
                posm = posm + first_aa_model
                posmf = posmf + first_aa_model

            discont = []
            for res in mdel_res_list:
                if is_aa(res) is True:
                    discont.append(res.get_id()[1])
                    last_aa_model = res.get_id()[1]
            a = 0
            for a in range(len(discont)-1):
                if discont[a] < posmf:
                    if discont[a+1] != discont[a] + 1:
                        posmf += 1
            # Define the ranges for the superimposer
            sup = Superimposer()
            ata_q1 = range(posq, posq + 40)
            ata_q2 = range(posqf, posqf + 40)
            ata_m1 = range(posm, posm + 40)
            ata_m2 = range(posmf, posmf + 40)
            # Create the list of the residues to be superimposed
            ref_atoms = []
            q_atoms = []
            # Add to the created lists the residues identified before
            try:
                for ref_chain in chainmodel[0]:
                    for ref_res in ref_chain:
                        if ref_res.get_id()[1] in ata_m1 or ref_res.get_id()[1] in ata_m2:
                            ref_atoms.append(ref_res['CA'])
                for q_chain in chainquery[0]:
                    for q_res in q_chain:
                        if q_res.get_id()[1] in ata_q1 or q_res.get_id()[1] in ata_q2:
                            q_atoms.append(q_res['CA'])
                    # Use the created lists to execute superimposer if they have the same number of atoms
                if len(ref_atoms) == len(q_atoms):
                    sup.set_atoms(ref_atoms, q_atoms)
                    sup.apply(chainquery)
                    io = PDBIO()
                    io.set_structure(chainquery)        
                    # Refine the model until rms>5 or 20 iterations
                    print("\nRefining the model for chain " + str(toalign[k_chain]) + "...")
                    time.sleep(2)
                    if numpy.abs(sup.rms) > 7.5:
                        checks[str(k_chain)] = 1
                        print("Careful, the superposition of the chains has an RMS of >5 (RMS = "
                              + str((numpy.abs(sup.rms)))[0:4] + ").\n\n")
                        io.save("query_aligned" + str(k_chain) + ".pdb")
                        time.sleep(1)
                    else:
                        print("Perfect for chain " + toalign[k_chain] + ".\n")
                        io.save("query_aligned" + str(k_chain) + ".pdb")
                        time.sleep(2)
                else:
                    print("\nThere is some problem with chain " + str(toalign[k_chain])
                          + ". Trying an alternate solution...")
                    time.sleep(2)
                    # Identify the most complete chain
                    chain_lengths = []
                    for chain in ordered_structure.get_chains():
                        chain_lengths.append(len(chain))
                    complete_chain = toalign[chain_lengths.index(max(chain_lengths))]
                    good_model = parser.get_structure("CH0" + complete_chain, "chain" + complete_chain + ".pdb")
                    try:
                        # Move a copy of the most complete chain to the position of the failed chain
                        range_model = range(int(last_aa_model/2), int((last_aa_model/2)+20))
                        fchain_atoms = []
                        goodm_atoms = []
                        for ref_chain in chainmodel[0]:
                            for ref_res in ref_chain:
                                if ref_res.get_id()[1] in range_model:
                                    fchain_atoms.append(ref_res['CA'])
                        for goodm_chain in good_model[0]:
                            for goodm_res in goodm_chain:
                                if goodm_res.get_id()[1] in range_model:
                                    goodm_atoms.append(goodm_res['CA'])      
                        sup.set_atoms(fchain_atoms, goodm_atoms)
                        sup.apply(good_model)
                        io = PDBIO()
                        io.set_structure(good_model)
                        io.save("good_model" + toalign[k_chain] + ".pdb")
                        # redefine the chainmodel to use the new superimposed file
                        chainmodel_fixed = parser.get_structure("CH0M", "good_model" + toalign[k_chain] + ".pdb")
                        ref_atomsfixed = []      
                        q_atoms = []
                        for m_chain in chainmodel_fixed[0]:
                            for m_res in m_chain:
                                if m_res.get_id()[1] in range_model:
                                    ref_atomsfixed.append(m_res['CA'])
                        for q_chain in chainquery[0]:
                            for q_res in q_chain:
                                if q_res.get_id()[1] in range_model:
                                    q_atoms.append(q_res['CA'])
                        # Superimpose query with the new structure
                        sup.set_atoms(ref_atomsfixed, q_atoms)
                        sup.apply(chainquery)
                        io.set_structure(chainquery)
                        print("Superimposing the query to a fixed chain " + str(toalign[k_chain]) + "...")
                        io.save("query_aligned" + str(k_chain) + ".pdb")
                        time.sleep(3)
                        count = 0
                        if numpy.abs(sup.rms) > 7.5:
                            checks[str(k_chain)] = 1
                            print("Careful, the superposition of the chains has an RMS of >7.5 (RMS = "
                                  + str((numpy.abs(sup.rms)))[0:4] + ").\n")
                            time.sleep(3)
                        else:
                            print("Perfect for chain " + toalign[k_chain] + ".")
                            time.sleep(3)
                    except BaseException as error:
                        checks[str(k_chain)] = 2
                        print("Didn't work either. Check the PDB file for discontinuties or errors.", error)
            except KeyError:
                checks[str(k_chain)] = 2
                print("Superimposition for chain " + str(toalign[k_chain]) + " did not work. "
                      "It could be that it's different from the initally modelled monomer.\n")
                break
        return checks

    chainid = []
    for chains in structure.get_chains():
        if chains.get_id() not in chainid:
            chainid.append(chains.get_id())

    number_of_chains = len(chainid)

    if number_of_chains > 1:
        i = 0
        for i in range(number_of_chains):
            shutil.copyfile("query_monomer.pdb", "query_aligned" + str(i) + ".pdb")
        # Give the newly created files a chain id.#
        chainid = []
        for chains in structure.get_chains():
            chainid.append(chains.get_id())
        chainrange = range(0, len(chainid))
        filenames = []
        for files in glob.glob("query_aligned*.pdb", recursive=True):
            filenames.append(files)
        i = 0
        for files in filenames:
            while i < len(filenames):
                mdl = model(env, file=filenames[i])
                mdl.rename_segments(segment_ids=[chainid[i]])
                mdl.write(file="query_aligned"+str(i)+".pdb")
                i += 1
        # split the chains avoiding any error if atoms disordered - bug in biopython
        # (https://www.biostars.org/p/380566/ and https://github.com/biopython/biopython/issues/455),
        # changed the following def in the PDB python file of biopython. This converst disordered residues into ordered
        # ones by eliminating all locations apart from the first one.#
        #############################################################################################
        # def get_unpacked_list(self):                                                               #
        #     """                                                                                   #
        #     Returns all atoms from the residue,                                                   #
        #     in case of disordered, keep only first alt loc and remove the alt-loc tag             #
        #     """                                                                                   #
        #     atom_list = self.get_list()                                                           #
        #     undisordered_atom_list = []                                                           #
        #     for atom in atom_list:                                                                #
        #         if atom.is_disordered():                                                          #
        #         if atom.is_disordered():                                                          #
        #             atom.altloc=" "                                                               #
        #             undisordered_atom_list.append(atom)                                           #
        #         else:                                                                             #
        #             undisordered_atom_list.append(atom)                                           #
        #     return undisordered_atom_list                                                         #
        #############################################################################################
        io = PDBIO()
        io.set_structure(structure)
        io.save("ord_pdb" + id_1.casefold() + ".ent")
        ordered_pdbfile = "ord_pdb" + id_1.casefold() + ".ent"
        file = id_1.casefold()
        ordered_structure = parser.get_structure(file, ordered_pdbfile)
        k = 0
        for chain in ordered_structure.get_chains():
            for k in range(number_of_chains):
                chain = ordered_structure[0][chainid[k]]
                io.set_structure(chain)
                io.save("chain" + str(chainid[k]) + ".pdb")
                k += 1

            # Use the method defined to superimpose the chains.
        toalign = {}
        ch = 0
        for chains in chainid:
            if ch < number_of_chains:
                toalign[ch] = chainid[ch]
                ch += 1
        checks = superimpose_chains(toalign)
        # Add the files together to create a query_template with all subunits
        filenames = glob.glob("query_aligned*.pdb")
        if 2 not in checks.values():
            if os.path.isfile("query_multimeric") == True:
                with open('quat_multimeric.pdb', 'w+') as fullquerymodel:  
                    for names in filenames: 
                        with open(names) as infile: 
                            for line in infile:
                                if 'END' not in line:
                                    fullquerymodel.write(line)
                    fullquerymodel.write("END\n")
            else:
                with open('quat_as_' + answer_pdb_id + '.pdb', 'w+') as fullquerymodel:  
                    for names in filenames: 
                        with open(names) as infile: 
                            for line in infile:
                                if 'END' not in line:
                                    fullquerymodel.write(line)
                    fullquerymodel.write("END\n")
            if 1 in checks.values():
                print("While creating your multimeric model some errors occured. You can find your model in "
                      + str(os.getcwd()) + " in a file named query_multimeric.pdb. BUT CHECK IT BEFORE USING IT IN ANY "
                                           "FURTHER STEP!")
            else:
                print("\tFinished running the module Blast&Modeller module!\nCreation of your model was succesfull."
                      "You can find it in " + str(os.getcwd()) + " in a file named query_multimeric.pdb along with "
                      "the monomeric modelled structure (query_monomeric.pdb), the sequences in fasta format, "
                      "the blast result and the alignment. You can continue now running the ActiveSiteID!")
        else:
            print("No multimeric assembly could be modelled due to errors in the chain or too high RMS. "
                  "Please check the pdb file for discontinuities or related problems "
                  "(check https://www.rcsb.org/structure/" + str(id_1) + "for more information). "
                  "You could try to rerun the stand alone multimeric builder with another PDB id.")
    else:
        mdl = model(env, file="query.B99990001.pdb")
        mdl.rename_segments(segment_ids="A")
        mdl.write(file="query.B99990001.pdb")
        print("\tFinished running the Quaternary submodule!\n You can find your result files in "
              + str(working_directory) + " in the Blast&Modeller folder.")

    # Clean the directory of intermediate files which are not necessary
    os.rmdir("obsolete")
    try:
        os.rename(pdbfile, 'quat_template.pdb')
    except BaseException as error:
        os.remove('quat_template.pdb')
        os.rename(pdbfile, 'quat_template.pdb')

    keep = ["quaternary_assembly.pdb", "query_multimeric.pdb", "query_monomer.pdb", "blast_result.xml", "quat_template.pdb", "query.fasta", "model.fasta", 'quat_as_' + answer_pdb_id + '.pdb']
    files = []
    for f in os.listdir("."):
        if os.path.isfile:
            files.append(f)

    for f in files:
        if f not in keep:
            os.remove(f)

    os.chdir(initial_location)
    

if __name__ == "__main__":
    main(folder_name, monomeric_model, multimeric_model)
