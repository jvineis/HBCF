#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Extract the cds of a genewise output and the allelic cds from the A.vaga genome and write both sequences to a single fasta file for alignment with muscle''')
parser.add_argument('fasta_file', help = 'the fasta file containing the sequence hit by the unpaired gene')
parser.add_argument('--genewise', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/GENEWISE/', help='The directory containing the genewise output.  Files within the directory must have names that link it to the fasta file input name e.g. fasta_file.genewise.out.txt')
parser.add_argument('--peptides', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/CDS_ALIGNMENT/', help='The location for output of the fasta files')
args = parser.parse_args()

genewise = args.fasta_file # name of the scaffold fasta file containing sequence hit by unpaired genes
genewise_in = open(args.genewise+genewise+'.genewise.out.txt', 'rU') # open the genewise output
outfile = open(args.peptides+genewise+'.cds_to_align.fa', 'w') # open the outfile to write the cds sequences
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


