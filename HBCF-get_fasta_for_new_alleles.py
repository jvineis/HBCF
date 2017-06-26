#!/bin/bash/env python
# This script takes the output from the HBCF-filter_blast_for_pair_correcting.py or any table with four columns where the 
# column ids are [unpaired_gene, unpaired_scaffold, new_gene_name, new_gene_scaffold, new_gene_start, new_gene_stop]
# The output are indiviudal fasta files for each row in the table with a name corresponding to the columns of the table.

from Bio import SeqIO
import sys


fasta = open('/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.scaffolds.fa', 'rU')
matching_locations = open(sys.argv[1], 'rU')# might be called NEW_GENES.txt contains columns "UNPAIRED_GENE, Unpaired scaffold, new_gene_name, new_gene_scaffold, new_gene_scaffold_start, new_gene_scaffold_end:

matching_list_fwd = []
matching_list_rev = []
with matching_locations as f:
    first_line = f.readline()
    for line in f:
        x = line.strip().split("\t")
        if int(x[4]) < int(x[5]):
            matching_list_fwd.append([x[0],x[1],x[2],x[3],x[4],x[5]])
        elif int(x[4]) > int(x[5]):
            matching_list_rev.append([x[0],x[1],x[2],x[3],x[5],x[4]])

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
            if int(scaf_start) > 201 and len(seq_dict[key]) > int(scaf_end)+200: # if the start of the hit is greater than 200nt from the start of the gene and the end is less than 200nt from the end  
                seq_range = (seq_dict[key][int(scaf_start)-200:int(scaf_end)+200])#add 200nt to the beginning and end of the sequence
                print("this seq is normal")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")#write the sequence to the fasta file
            elif int(scaf_start) < 201 and len(seq_dict[key]) > int(scaf_end)+200:# if the start of the hit is within the first 200bp of the scaffold seq and the end is prior to the last 200bp of the scaffold
                seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)+200])
                print("couldn't write the 200bp front end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
            elif int(scaf_start) > 201 and len(seq_dict[key]) < int(scaf_end)+200:# if the end of the hit is within the last 200bp of the scaffold seq and start is beyond the first 200bp of the scaffold
                seq_range = (seq_dict[key][int(scaf_start)-200:int(scaf_end)])
                print("couldn't write the 200bp tail end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
            elif int(scaf_start) < 201 and len(seq_dict[key]) < int(scaf_end)+200:# if the end and beginning of the hit are both outside the 200bp buffer that is being added to the sequence
                seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)])
                print("couldn't write the 200bp tail or end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
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
            if int(scaf_start) > 201 and len(seq_dict[key]) > int(scaf_end)+200: # if the start of the hit is greater than 200nt from the start of the gene and the end is less than 200nt from the end
                seq_range = (seq_dict[key][int(scaf_start)-200:int(scaf_end)+200].reverse_complement())#add 200nt to the beginning and end of the sequence
                print("rev-this seq is normal")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")#write the sequence to the fasta file
            elif int(scaf_start) < 201 and len(seq_dict[key]) > int(scaf_end)+200:# if the start of the hit is within the first 200bp of the scaffold seq and the end is prior to the last 200bp of the scaffold
                seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)+200].reverse_complement())
                print("rev-couldn't write the 200bp front end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
            elif int(scaf_start) > 201 and len(seq_dict[key]) < int(scaf_end)+200:# if the end of the hit is within the last 200bp of the scaffold seq and start is beyond the first 200bp of the scaffold
                seq_range = (seq_dict[key][int(scaf_start)-200:int(scaf_end)].reverse_complement())
                print("rev-couldn't write the 200bp tail end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
            elif int(scaf_start) < 201 and len(seq_dict[key]) < int(scaf_end)+200:# if the end and beginning of the hit are both outside the 200bp buffer that is being added to the sequence
                seq_range = (seq_dict[key][int(scaf_start):int(scaf_end)].reverse_complement())
                print("rev-couldn't write the 200bp tail or end buffer")
                outfile.write(">"+str(outfile_header)+"\n"+ str(seq_range) + "\n")
