#!/bin/bash/env python

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(description = '''This script takes the output from the HBCF-filter_blast_for_pair_correcting.py or any table with four columns where the column ids are [unpaired_gene, unpaired_scaffold, new_gene_name, new_gene_scaffold, new_gene_start, new_gene_stop]. The output are indiviudal fasta files for each row in the table with a name corresponding to the columns of the table.  If the start position is larger than the stop position, the sequence returned is the reverse complement''')
parser.add_argument('--fasta_file', default = '/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.scaffolds.fa', help = 'a fasta file to extract sequences from')
parser.add_argument('matching_locations', help = 'a table with the [unpaired_gene, unpaired_scaffold, new_gene_name, new_gene_scaffold, new_gene_start, new_gene_stop] columns')
args = parser.parse_args()

fasta = open(args.fasta_file, 'rU')
matching_locations = open(args.matching_locations, 'rU')# contains columns Upaired_gene, Unpaired scaffold, new_gene_name, new_gene_scaffold, new_gene_scaffold_start, new_gene_scaffold_end:

matching_list_fwd = []
matching_list_rev = []
with matching_locations as f:
    first_line = f.readline()
    for line in f:
        x = line.strip().split("\t")
        if int(x[4]) < int(x[5]):
            matching_list_fwd.append([x[0],x[1],x[2],x[3],x[4],x[5]])# this ensuers that the start and stop are correct
        elif int(x[4]) > int(x[5]):
            matching_list_rev.append([x[0],x[1],x[2],x[3],x[5],x[4]])# this ensures that the start and stop are correct

seq_dict = {}
for record in SeqIO.parse(fasta, "fasta"):
    #header = record.id.split(" ")
    seq_dict[record.id] = record.seq

for line in matching_list_fwd:
    scaf_name =  line[3]
    scaf_start = line[4]
    scaf_end = line[5]
    new_name = line[2]
    outfile_header = new_name+":"+scaf_name+":"+scaf_start+":"+scaf_end
    outfile_name = new_name+"_"+scaf_name+"_"+scaf_start+"_"+scaf_end+".fa"
    print (outfile_header)
    for key in seq_dict.keys():
        header = key.split(" ")
        if scaf_name == header[0]:
            outfile = open(outfile_name, 'w')
            seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)])#defined by the start and stop positions of the blast hit
            outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")#write the sequence to the fasta file

for line in matching_list_rev:
    scaf_name =  line[3]
    scaf_start = line[4]
    scaf_end = line[5]
    new_name = line[2]
    outfile_header = new_name+":"+scaf_name+":"+scaf_start+":"+scaf_end
    outfile_name = new_name+"_"+scaf_name+"_"+scaf_start+"_"+scaf_end+".fa"
    print (outfile_header)
    for key in seq_dict.keys():
        header = key.split(" ")
        if scaf_name == header[0]:
            outfile = open(outfile_name, 'w')
            seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)])#defined by the start and stop positions of the blast hit
            outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")#write the sequence to the fasta file
