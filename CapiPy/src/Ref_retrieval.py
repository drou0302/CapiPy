#!/usr/bin/env python3

# Immobilization reference retrieval v.0.1
# Copyleft David Roura, 2020
"""Retrieval of 20 most relevant references using the keywords based on the UniProt entry of the query protein

"""
import time
import urllib.request
import urllib.parse
import os
import shutil
import csv

from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio import SeqIO, SearchIO
from metapub import PubMedFetcher

# Locate the workspace#


def main(folder_name, query):
    initial_location = os.getcwd()
    working_directory = folder_name
    if os.path.exists(working_directory) is True:
        os.chdir(working_directory)
        # Check if ActiveSite directory exists
        if os.path.exists("RelevantPapers") is False:
            os.mkdir("RelevantPapers")
        os.chdir("./RelevantPapers")
    else:
        os.mkdir(working_directory)
        os.chdir(working_directory)
        os.mkdir("RelevantPapers")
        os.chdir("./RelevantPapers")


# Check if the files from the Blast&Modeller exist
    if query == 'YES':    
        if os.path.exists("../Blast&Modeller") is True:
            shutil.copy("../Blast&Modeller/query.fasta", "./query.fasta")
        else:
            print("The query file doesn't exists! Make sure you have run the BLAST&Modeller and the query.fasta is ok!")
    else:
        with open("query.fasta", "w+") as f:
            protein_sequence = query
            f.write(">query" + "\n" + protein_sequence)

# Using the query sequence copied or indicated by the user, perform a blast to identify the uniprot identifier
    try:
        query = SeqIO.read(open("query.fasta"), format="fasta")
        print("Blast search running online... This migth take a while.")
        result_handle = NCBIWWW.qblast("blastp", "swissprot", query, auto_format="XML", matrix_name="BLOSUM62",
                                       expect=0.0001, word_size="6", gapcosts="11 1", alignments=10)
        print("Blast succesfull!")
        with open("blast_result.xml", "w+") as blast_result:
            blast_result.write(result_handle.read())
    except BaseException as ex: 
        print("Some error occured:" + ex.message)
        time.sleep(5)
        quit()

# Check if blast was succesfull
    try:
        blastup = SearchIO.read("blast_result.xml", "blast-xml")
    except BaseException as ex:
        print("Some error occured during your blast search.\n" + ex.message)
        time.sleep(5)
        quit()

# Extract uniprotID and the protein name
    uniprot_id = blastup[0].id.split("|")[1][0:6]
    if uniprot_id.isalpha():
        uniprot_id = input("It seems there is some problem finding the Uniprot ID of your protein. "
                           "If you have it, please enter it. Else, press enter to exit.")
    else:
        print("Uniprot id found (" + str(uniprot_id) + "). Extracting information...")
    handle = urllib.request.urlopen("https://www.uniprot.org/uniprot/"+str(uniprot_id)+".xml")
    record = SeqIO.read(handle, "uniprot-xml")
    with open("uniprot_papers.txt", "w+") as papers1:
        for papers in record.annotations["references"]:
            papers1.write(str(papers) + "\n")

    protein = []

    for info in record.annotations:
        if info == "submittedName_fullName":
            protein.append(record.annotations[info])
        elif info == "recommendedName_fullName":
            protein.append(record.annotations[info])
        elif info == "alternativeName_fullName":
            protein.append(record.annotations[info])

# Use the protein name, term by term, as keywords and add immob*
    keywords = []
    flat_keywords = []
    for names in protein:
        for name in names:
            keywords.append(name.split(" "))
    
    for list_1 in keywords:
        for a in list_1:
            flat_keywords.append(a)
        
    keywords = ""
    for words in flat_keywords:
        keywords += (str(words) + " OR ")

    keywords = keywords[0:-3]
    keywords += "AND immob*"

    # Search on pubmed using Entrez from Biopython
    print("Search on PubMed database is going to start with: \"" + keywords + "\" as the keywords.")
    Entrez.email = "example@example.com"
    handle1 = Entrez.esearch(db="pubmed", 
                             sort="relevance",
                             retmax="20",
                             retmode="xml",
                             term=keywords)
    article_list = Entrez.read(handle1)
    article_id_list = article_list["IdList"]

    uniprot_articles = []
    with open("uniprot_papers.txt", "r") as papers1:
        for line in papers1:
            if "pubmed id" in line:
                uniprot_articles.append(line[11:-1])
    for i in range(len(uniprot_articles)-1):
        try:
            if uniprot_articles[i] == "":
                del(uniprot_articles[i])
        except ErrorBase:
            break

    # Use metapub package to retreive the information and write it in a csv file
    print("Retrieving the information...")
    with open("relevant_papers.csv", "w", newline="", encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(["Number", "Article ID", "Title", "Year", "Link", "DOI"])
        fetcher = PubMedFetcher()
        for i in range(0, len(article_id_list)-1):
            src = fetcher.article_by_pmid(article_id_list[i])
            number = i+1
            article_id = article_id_list[i]
            title = src.title
            year = src.year
            link = "https://pubmed.ncbi.nlm.nih.gov/" + article_id_list[i]
            DOI = src.doi
            writer.writerow([number, article_id, title, year, link, DOI])
        for i in range(0, len(uniprot_articles)-1):
            up_src = fetcher.article_by_pmid(uniprot_articles[i])
            number = "Uniprot" + str(i+1)
            article_id = uniprot_articles[i]
            title = up_src.title
            year = up_src.year
            link = "https://pubmed.ncbi.nlm.nih.gov/" + uniprot_articles[i]
            DOI = up_src.doi
            writer.writerow([number, article_id, title, year, link, DOI])

    print("\tFinished running the module Reference retrieval module!\n You can find your result files in "
          + str(working_directory) + " in the RelevantPapers folder.\n The papers are organized in the csv file named "
                                     "\"relevant_papers.csv!\"")
    os.chdir(initial_location)

if __name__ == "__main__":
    main(folder_name, query)
