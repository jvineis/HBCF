#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio import pairwise2
import csv
import StringIO

parser = argparse.ArgumentParser(description='''This script creates two fasta files representing the matching allele pairs of the Adineta vaga genome.  They can be useful for mapping in order to detect SNPs that exist between pairs''')
parser.add_argument('-pairs_txt_file', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles', help = 'The allele pairs txt file containing the genes that are identified as alleles.  Only the 3rd and 4th column containing the gene names eg [GSADVT00000227001]    [GSADVT00022299001] is used from this file.')
parser.add_argument('-fasta_file', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa', help = 'A fasta file that contains gene names in the header eg. [>GSADVT00000001001]')
args = parser.parse_args()


handle = open(args.fasta_file, 'rU')
pairs = open(args.pairs_txt_file, 'rU')


pair_dict = {}
seq_list1 = []
seq_list2 = []
for line in pairs:
    x = line.rstrip("\n").split("\t")
    seq_list1.append(x[0])
    seq_list2.append(x[1])


pair_list1 = []
pair_list2 = []

for rec in SeqIO.parse(handle, "fasta"):
    if rec.id in seq_list1:
        pair_list1.append(rec)
    if rec.id in seq_list2:
        pair_list2.append(rec)
output_handle1 = open("allele_pairs1.fa", "w")
output_handle2 = open("allele_pairs2.fa", "w")
SeqIO.write(pair_list1, output_handle1, "fasta")
SeqIO.write(pair_list2, output_handle2, "fasta")
output_handle1.close()
output_handle1.close()


