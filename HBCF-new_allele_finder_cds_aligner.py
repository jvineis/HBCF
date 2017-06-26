#!/usr/bin/env python

import sys
from Bio import SeqIO

# This script is designed to align the CDS of trancsripts predicted by genewise with the CDS found in the 
# /groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa fasta 
# The script can be run from anywhere in the server and will put the output in the "CDS_ALIGNMENT" dicectory
# indicated below


genewise = sys.argv[1] # name of the scaffold fasta file containing sequence hit by unpaired genes
genewise_in = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/GENEWISE/'+genewise+'.genewise.out.txt', 'rU') # open the genewise output
outfile = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/CDS_ALIGNMENT/'+genewise+".cds_to_align.fa", 'w') # open the outfile to write the cds sequences
file_name = genewise.split("_") # separating the characters of the fasta file neede to select from the fasta containg the cds
print file_name
gsadv = file_name[0][0:14]+'001' # the characters needed to select the cds from the Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa
print(gsadv)
cds_fasta = open('/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa', 'rU')

with cds_fasta as f:
    cds = SeqIO.parse(f, "fasta")
    for rec in cds:
        if gsadv in rec.id:
            outfile.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n')

bases = ['A','G','C','T']

sequence = []
header = []
for line in genewise_in:
    x = line.strip()
    if "sp" in x:
        header.append(x)
    if x[0] in bases:
        sequence.append(x)
        
outfile.write(str(header[0])+'\n')

for seq in sequence:
    outfile.write(str(seq))


