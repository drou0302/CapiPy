#!/usr/bin/env python3

"""#Blast&modeller v.0.1
'#Copyleft David Roura, 2020
Combination of BLAST and MODELLER softwares with Superimposer to create models of monomeric and multimeric proteins

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


def main(folder_answer, folder_name, prot_sequence, blast_response, clustalw_exe):
    initial_location = os.getcwd()
    """Creation of the working directory"""
    work_creation = folder_answer
    if work_creation == "NO":
        working_directory = folder_name
        if os.path.exists('./' + working_directory) is True:
            os.chdir(working_directory)
            if os.path.exists("Blast&Modeller") is False:
                os.mkdir("Blast&Modeller")
                os.chdir("./Blast&Modeller")
            else:
                os.chdir("./Blast&Modeller")
        else:
            print("Error accessing the folder. Are you sure it exists?")
            quit()
            
    elif work_creation == "YES":
        working_directory = folder_name
        if os.path.exists("./" + working_directory) is False:
            os.mkdir("./" + working_directory)
            print("Directory ", working_directory, " created.")
            os.chdir("./" + working_directory)
        else:
            os.chdir("./" + working_directory)
        os.mkdir("Blast&Modeller")
        os.chdir("./Blast&Modeller")

    # Input of the protein sequence

    protein_sequence = prot_sequence

    with open("query.fasta", "w+") as f:
        f.write(">query" + "\n" + protein_sequence)

    # Eliminate the his tag from the search to avoid partial hits with no real significance#

    if "HHHHH" in protein_sequence[0:20]:
        with open("protein_blast.fasta", "w+") as f:
            f.write(">query\n" + protein_sequence[20:])
    elif "HHHHH" in protein_sequence[len(protein_sequence) - 20:len(protein_sequence)]:
        with open("protein_blast.fasta", "w+") as f:
            f.write(">query\n" + protein_sequence[0:len(protein_sequence) - 20])
    else:
        with open("protein_blast.fasta", "w+") as f:
            f.write(">query" + "\n" + protein_sequence)

    # Identification of the best hit from the PDB database using query.fasta as query in blast
    blast_input = blast_response
    if blast_input == "WEB":
        query = SeqIO.read(open("protein_blast.fasta"), format="fasta")
        print("Blast search running online...")
        result_handle = NCBIWWW.qblast("blastp", "pdb", query, auto_format="XML", matrix_name="BLOSUM62",
                                       expect=0.0001, word_size="6", gapcosts="11 1", alignments=10)
        print("Blast succesfull!")
        with open("blast_result.xml", "w+") as blast_result:
            blast_result.write(result_handle.read())
    elif blast_input == "LOCAL":
        result_handle = NcbiblastpCommandline(query="protein_blast.fasta", db="pdbaa", evalue=0.0001, outfmt=5,
                                              out="blast_result.xml")
        result_handle()
        print("Blast search running locally...")

        # Extract first hit and create a fasta with its pdb identifier and its sequence

    blast_qresult = SearchIO.read("blast_result.xml", "blast-xml")
    if len(blast_qresult) > 0:
        blast_qresult = SearchIO.read("blast_result.xml", "blast-xml")[0]
        best_hit = blast_qresult[0]
        id_1 = best_hit.hit.id.split("|")[1]
        print("Your blast returned the following as the best hit:\n" + str(best_hit)
              + ".\n Visit https://www.rcsb.org/structure/" + str(id_1) + " for more information about the model.")
        time.sleep(3)
    else:
        print("Your blast doesn't seem to contain any result. Check the sequence and run the program again. "
              "If you run it over the server, consider using local blast. The program is going to quit now.")
        time.sleep(3)
        quit()

    # Download from the PDB the model

    pdbl = PDBList()
    pdbl.retrieve_pdb_file(id_1, pdir=".", file_format="pdb")
    parser = PDBParser()
    ppb = PPBuilder()
    pdbfile = "pdb" + id_1.casefold() + ".ent"
    file_name = str(id_1.casefold())
    structure = parser.get_structure(file_name, pdbfile)
    chain = structure.get_chains()
    env = environ()

    # Method for alignment with modeller. The option determines if only aa are used for the alignment or not.
    # This helps for enzymes with discontinuities and internal cofactors such as MIO dependant enzymes.
    def align_for_modeller(option):
        for chain in structure.get_chains():
            last_chain = chain.id
            tnum_res = (len([_ for _ in chain.get_residues() if PDB.is_aa(_)]))
        for pp in ppb.build_peptides(structure[0][str(last_chain)], aa_only=option):
            with open("model.fasta", "a") as f:
                sequence = pp.get_sequence()
                f.write(str(sequence))
        mdel = structure[0]
        chain = mdel[last_chain]
        res_list = PDB.Selection.unfold_entities(chain, "R")
        for residue in res_list:
            nid_1aa = (res_list[0].get_id())[1]
            nid_lastaa = (res_list[tnum_res - 1].get_id())[1]
        # Modify files so they have only 1 identifier + seq (id_1 or query)#
        with open('model.fasta', 'r') as o:
            data = o.read()
        with open('model.fasta', 'w') as mo:
            mo.write(">" + file_name + "\n" + data)
        with open("query.fasta") as o2:
            line2 = o2.readlines()
        line2[0] = ">query\n"
        with open("query.fasta", "w") as m1:
            m1.writelines(line2)
            # Combine original sequence and first hit into a fasta input file for the alignment#
        filenames = ["query.fasta", "model.fasta"]
        with open('alignment.fasta', 'w+') as aligninput:
            for files in filenames:
                with open(files) as infile:
                    aligninput.write(infile.read())
                aligninput.write("\n")
                # Profile-profile alignment using salign from modeller#
        log.none()
        aln = alignment(env, file='alignment.fasta', alignment_format='FASTA')
        aln.salign(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, 0), output='', align_block=15,
                   align_what='PROFILE', alignment_type='PAIRWISE', comparison_type='PSSM', similarity_flag=True,
                   substitution=True, smooth_prof_weight=10.0)
        aln.write(file='salign.ali', alignment_format='PIR')
        print("Alignment of template and query for modelling successfull.")
        # Fix formatting of the .ali file to specify the 1st and last aminoacid and the structure of the model protein#
        shutil.copyfile("salign.ali", "salign1.ali")
        with open("salign1.ali", "r") as file:
            filedata = file.read()
        replacement = filedata.replace(">P1;" + str(file_name) + "\nsequence::     : :     : :::-1.00:-1.00",
                                       ">P1;" + str(file_name) + "\nstructureX" + ":" + str(file_name) + ":" + str(
                                           nid_1aa) + ":" + str(last_chain) + ":" + str(nid_lastaa) + ": ::::")
        with open("salign1.ali", "w+") as f:
            f.write(replacement)
        # Make a single model of the query sequence using: salign1.ali and model.pdb#
        a = automodel(env, alnfile="salign1.ali", knowns=id_1.casefold(), sequence="query")
        a.starting_model = 1
        a.ending_model = 1
        a.make()
        # Check how good the alignment is#
        pir_alignment = AlignIO.read("salign1.ali", "pir")
        total_length = len(pir_alignment[0])
        gaps_1 = 0
        gaps_2 = 0
        for aas in pir_alignment[0].seq:
            if aas == "-":
                gaps_1 += 1
        for aas in pir_alignment[1].seq:
            if aas == "-":
                gaps_2 += 1
        if gaps_1 / total_length > 0.5 or gaps_2 / total_length > 0.5:
            print(
                "\nYour model protein covers less than half of the query. The created model could be inaccurate, " +
                "please check the result before continuing with anything else.\n")
            time.sleep(3)

    # Run the previous method, changing the option if it fails the first time. If non of them work,
    # raise and error and quit the program.
    try:
        align_for_modeller(True)
        print("First model created succesfully")
        time.sleep(3)
    except BaseException:
        try:
            delete_files = ["alignment.fasta", "salign.ali", "salign1.ali", "model.fasta"]
            for files in delete_files:
                if os.path.exists(files):
                    os.remove(files)
            align_for_modeller(False)
        except BaseException:
            delete_files = ["salign.ali", "salign1.ali", "model.fasta"]
            for files in delete_files:
                if os.path.exists(files):
                    os.remove(files)
            print(
                "\nThere is problem with your model. It could be related to unnatural amino acids or problems in "
                "the residue numbering. Please, retry with a different PDB or try and fix the errors in the file."
                "\nThe program is now going to quit. Sorry!")
            quit()

    # Method to superimpose the chains#
    def superimpose_chains(dict):
        checks = {}
        # Check if clustalw is available.
        for k_chain in toalign.keys():
            chainmodel = parser.get_structure("CH0" + toalign[k_chain], "chain" + toalign[k_chain] + ".pdb")
            chainquery = parser.get_structure("Q00" + str(k_chain), "query_aligned" + str(k_chain) + ".pdb")
            ppb = PPBuilder()
            # To maximise the accuracy of the superimposer and avoid problems with alignments of low quality instead
            # of using the whole protein to superimpose, two regions at the start and end will be used. To do so, first
            # it's necessary to identify two regions at the start and end where both there are no gaps, checking the
            # alignment for a minimum of 40 positions without "-" in reading frames of 5 aminoacids.
            frac_chainmodelseq = []
            chainmodelseq = ""
            for pp in ppb.build_peptides(chainmodel):
                frac_chainmodelseq.append(pp.get_sequence())
            for fractions in frac_chainmodelseq:
                chainmodelseq += fractions
            for pp in ppb.build_peptides(chainquery):
                chainqueryseq = pp.get_sequence()
            with open("initial" + str(k_chain) + ".fasta", "w") as input_file:
                input_file.write(
                    ">model_" + toalign[k_chain] + "\n" + str(chainmodelseq) + "\n>query_" + str(k_chain) + "\n" + str(
                        chainqueryseq))

            ex_clustal = ClustalwCommandline(clustalw_exe, infile="initial" + str(k_chain) + ".fasta")
            stdout, stderr = ex_clustal()
            alignment = AlignIO.read("initial" + str(k_chain) + ".aln", "clustal")
            length = len(alignment[0])
            # Identify the first region of 40 aa without gaps in the alignment
            i = 0
            posm = 0
            posq = 0
            posmf = 0
            posqf = 0
            for a in range(int(length / 5)):
                if "-" not in alignment[0][i:i + 40].seq and "-" not in alignment[1][i:i + 40]:
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
            k = length - (length // 4)
            while k < length - 40:
                if "-" not in alignment[0][k:k + 40].seq and "-" not in alignment[1][k:k + 40].seq:
                    posqf = k
                    posmf = k
                    break
                else:
                    k += 5
            if k + 40 > length:
                k = length - (length // 4)
                while k < length - 40:
                    if "-" not in alignment[0][k:k + 40].seq and "-" not in alignment[1][k:k + 40].seq:
                        posqf = k
                        posmf = k
                        break
                    else:
                        k += 2
            if posqf == 0:
                posqf = length - (length // 4)
                posmf = length - (length // 4)
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
            for a in range(len(discont) - 1):
                if discont[a] < posmf:
                    if discont[a + 1] != discont[a] + 1:
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
                        print("Careful, the superposition of the chains has an RMS of >5 (RMS = " + str(
                            (numpy.abs(sup.rms)))[0:4] + ").\n\n")
                        io.save("query_aligned" + str(k_chain) + ".pdb")
                        time.sleep(1)
                    else:
                        print("Perfect for chain " + toalign[k_chain] + ".\n")
                        io.save("query_aligned" + str(k_chain) + ".pdb")
                        time.sleep(2)
                else:
                    print("\nThere is some problem with chain " + str(
                        toalign[k_chain]) + ". Trying an alternate solution...")
                    time.sleep(2)
                    # Identify the most complete chain
                    chain_lengths = []
                    for chain in ordered_structure.get_chains():
                        chain_lengths.append(len(chain))
                    complete_chain = toalign[chain_lengths.index(max(chain_lengths))]
                    good_model = parser.get_structure("CH0" + complete_chain, "chain" + complete_chain + ".pdb")
                    try:
                        # Move a copy of the most complete chain to the position of the failed chain
                        range_model = range(int(last_aa_model / 2), int((last_aa_model / 2) + 20))
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
                            print("Careful, the superposition of the chains has an RMS of >7.5 (RMS = " + str(
                                (numpy.abs(sup.rms)))[0:4] + ").\n")
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

    # Check if the model protein is formed by  more than one chain. Although not prove of it being multimeric,
    # better more work than overwatch this. If so, ask user if a multimeric model using the template has to be created#
    chainid = []
    for chains in structure.get_chains():
        if chains.get_id() not in chainid:
            chainid.append(chains.get_id())

    number_of_chains = len(chainid)

    if number_of_chains > 1:
        i = 0
        for i in range(number_of_chains):
            shutil.copyfile("query.B99990001.pdb", "query_aligned" + str(i) + ".pdb")

        # Give the newly created files a chain id.#
        chainid = []
        for chains in structure.get_chains():
            chainid.append(chains.get_id())
        chainrange = range(0, len(chainid))
        filenames = []
        for files in glob.glob("query_aligned*.pdb", recursive=True):
            if files != "query.B99990001.pdb":
                filenames.append(files)
        i = 0
        for files in filenames:
            while i < len(filenames):
                mdl = model(env, file=filenames[i])
                mdl.rename_segments(segment_ids=[chainid[i]])
                mdl.write(file="query_aligned" + str(i) + ".pdb")
                i += 1
        # split the chains avoiding any error if atoms disordered
        # - bug in biopython (https://www.biostars.org/p/380566/ and https://github.com/biopython/biopython/issues/455),
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
        k = 1
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
            if ch < number_of_chains - 1:
                toalign[ch] = chainid[ch]
                ch += 1
        checks = superimpose_chains(toalign)
        # Add the files together to create a query_template with all subunits
        filenames = glob.glob("query_aligned*.pdb")
        if 2 not in checks.values():
            with open('query_multimeric.pdb', 'w+') as fullquerymodel:
                for names in filenames:
                    with open(names) as infile:
                        for line in infile:
                            if 'END' not in line:
                                fullquerymodel.write(line)
                fullquerymodel.write("END\n")
            if 1 in checks.values():
                print("While creating your multimeric model some errors occured. You can find your model in " + str(
                    os.getcwd()) + " in a file named query_multimeric.pdb. "
                                   "BUT CHECK IT BEFORE USING IT IN ANY FURTHER STEP!")
            else:
                print(
                    "\tFinished running the module Blast&Modeller module!\nCreation of your model was succesfull. "
                    "You can find it in " + str(os.getcwd()) + " in a file named query_multimeric.pdb along with the "
                    "monomeric modelled structure (query_monomeric.pdb), the sequences in fasta format, the "
                    "blast result and the alignment. You can continue now running the ActiveSiteID!")
    else:
        mdl = model(env, file="query.B99990001.pdb")
        mdl.rename_segments(segment_ids="A")
        mdl.write(file="query.B99990001.pdb")
        print("\tFinished running the module Blast&Modeller module!\n You can find your result files in " + str(
            working_directory) + " in the Blast&Modeller folder.\nIf you think your protein should have a multimeric " +
            "assembly, check https://www.rcsb.org/structure/" + str(id_1) +
            " for more information. You can continue now running the ActiveSiteID!")

    # Clean the directory of intermediate files which are not necessary
    os.rmdir("obsolete")
    os.rename("query.B99990001.pdb", "query_monomer.pdb")
    os.rename("salign1.ali", "alignment.ali")
    os.rename(pdbfile, 'template.pdb')
    keep = ["query_multimeric.pdb", "query_monomer.pdb", "blast_result.xml", 'template.pdb', "query.fasta", "model.fasta"]
    files = []
    for f in os.listdir("."):
        if os.path.isfile:
            files.append(f)

    for f in files:
        if f not in keep:
            os.remove(f)

    os.chdir(initial_location)


if __name__ == "__main__":
    main(folder_answer, folder_name, prot_sequence, blast_response)
