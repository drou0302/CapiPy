# CapiPy
CapiPy (Computer Assistance for Protein Immobilisation – Python) is a small collection of 4 main packages to retrieve useful information for the immobilisation of a query protein.	 
Information from each module, specially the Blast and modeller is preferentially used by the others. In case no such information is given, the user will have to input the necessary information. 

MODULE 1: BLAST and MODELLER: this module only requires the 1 letter amino acid sequence of the query protein. This will be blasted against the protein data base to identify the protein with the highest similarity and a model will be created. If the protein is multimeric, the option to model the quaternary structure is also given via superposition of the monomeric model onto the other chains. 

MODULE 2: ACTIVE SITE ID: the identification of the active site from the query protein is performed using the information present in UNIPROT and/or the M-CSA database. If the protein is annotated to use any cofactor, this information is also given. 

MODULE 3: SIZE AND CLUSTERS IDENTIFICATION: this is the most important module regarding immobilisation. First, the surface and volume of the model of the query protein (or a user specified PDB) are calculated. Secondly, all protein residues are classified according to their solvent exposure and from the exposed ones, clusters formed by residues important for immobilisation identified. The output can easily be visualised in PyMOL.

MODULE 4: IMMOBILISATION PAPER RETRIEVAL: to complement the information given by the other modules, in this last module 20 immobilisation papers related to the query protein are retrieved. To do so, the query used in the previous modules (or input here by the user) is blasted against the SwissProt database to identify the best hit.  The information of this best hit, such as protein name or identifiers, is extracted and used as keywords to submit a search into the Pubmed database. The 20 most relevant publications (by citation) are listed in an excel file.  

In addition, two stand-alone functionalities of the first and third module are also available:

MODULE 1.1: QUATERNARY STRUCTURE DETERMINATION: with a monomeric model created, this module allows the user to select a template from the PDB dataset to create a different quaternary assembly.

MODULE 3.1: CLUSTER DISTANCE: once the clusters have been identified, similar to the last part of the third module, this allows the calculation of the distance between the clusters and any user-specified position in the query protein. 
