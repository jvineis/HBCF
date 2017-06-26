#!/bin/usr/env python

import sys
from Bio import SeqIO
from Bio import Seq

unpaired_coords = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/UNPAIRED_GENE_COORDINATES.txt', 'rU')
new_coords = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_GENES.txt', 'rU')
exist_coords = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/EXISTING_GENES.txt', 'rU')

fasta = open('/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.scaffolds.fa', 'rU')

unpaired_list = []
with unpaired_coords as f:
    firstline = f.readline()
    for line in f:
        x = line.strip().split("\t")
        if int(x[2]) < int(x[3]):
            unpaired_list.append(x)
        elif int(x[2]) > int(x[3]):
            unpaired_list.append([x[0],x[1],x[3],x[2]])

fwd_new_list = []
rev_new_list = []
with new_coords as f:
    firstline = f.readline()
    for line in new_coords:
        x = line.strip().split("\t")
        if int(x[4]) < int(x[5]):
            fwd_new_list.append(x)
        elif int(x[4]) > int(x[5]):
            rev_new_list.append([x[0],x[1],x[2],x[3],x[5],x[4]])

fwd_exist_list = []
rev_exist_list = []
with exist_coords as f:
    firstline = f.readline()
    for line in exist_coords:
        x = line.strip().split("\t")
        if int(x[4]) < int(x[5]):
            fwd_exist_list.append(x)
        elif int(x[4]) > int(x[5]):
            rev_exist_list.append([x[0],x[1],x[2],x[3],x[5],x[4]])
            

seq_dict = {}
for record in SeqIO.parse(fasta, "fasta"):
    header = record.id.split(" ")
    header = header[0]
    seq_dict[header] = record.seq

for hit in rev_new_list:
    outname = hit[0]+'_'+hit[1]+'_'+hit[2]+'_'+hit[3]+'_'+hit[4]+'_'+hit[5]+'.fa'
    outfile = open(outname, 'w')
    for gene in unpaired_list:
        if hit[0] == gene[0] and hit[1] == gene[1]:
            hit_seq = seq_dict[hit[3]][int(hit[4]):int(hit[5])]
            hit_seq = hit_seq.reverse_complement()
            hit_header = hit[0]+":"+hit[2]+":"+hit[3]+":"+hit[4]+":"+hit[5]
            gen_seq = seq_dict[gene[1]][int(gene[2]):int(gene[3])]
            gen_header = gene[0]+":"+gene[1]+":"+gene[2]+":"+gene[3]
            outfile.write(">"+str(hit_header)+"\n"+str(hit_seq)+"\n"+">"+str(gen_header)+"\n"+str(gen_seq))

for hit in rev_exist_list:
    outname = hit[0]+'_'+hit[1]+'_'+hit[2]+'_'+hit[3]+'_'+hit[4]+'_'+hit[5]+'.fa'
    outfile = open(outname, 'w')
    for gene in unpaired_list:
        if hit[0] == gene[0] and hit[1] == gene[1]:
            hit_seq = seq_dict[hit[3]][int(hit[4]):int(hit[5])]
            hit_seq = hit_seq.reverse_complement()
            hit_header = hit[0]+":"+hit[2]+":"+hit[3]+":"+hit[4]+":"+hit[5]
            gen_seq = seq_dict[gene[1]][int(gene[2]):int(gene[3])]
            gen_header = gene[0]+":"+gene[1]+":"+gene[2]+":"+gene[3]
            outfile.write(">"+str(hit_header)+"\n"+str(hit_seq)+"\n"+">"+str(gen_header)+"\n"+str(gen_seq))

for hit in fwd_new_list:
    outname = hit[0]+'_'+hit[1]+'_'+hit[2]+'_'+hit[3]+'_'+hit[4]+'_'+hit[5]+'.fa'
    outfile = open(outname, 'w')
    for gene in unpaired_list:
        if hit[0] == gene[0] and hit[1] == gene[1]:
            hit_seq = seq_dict[hit[3]][int(hit[4]):int(hit[5])]
            hit_header = hit[0]+":"+hit[2]+":"+hit[3]+":"+hit[4]+":"+hit[5]
            gen_seq = seq_dict[gene[1]][int(gene[2]):int(gene[3])]
            gen_header = gene[0]+":"+gene[1]+":"+gene[2]+":"+gene[3]
            outfile.write(">"+str(hit_header)+"\n"+str(hit_seq)+"\n"+">"+str(gen_header)+"\n"+str(gen_seq))
            print(hit_header)
            print(hit_seq)
            print(gen_header)
            print(gen_seq)

for hit in fwd_exist_list:
    outname = hit[0]+'_'+hit[1]+'_'+hit[2]+'_'+hit[3]+'_'+hit[4]+'_'+hit[5]+'.fa'
    outfile = open(outname, 'w')
    for gene in unpaired_list:
        if hit[0] == gene[0] and hit[1] == gene[1]:
            hit_seq = seq_dict[hit[3]][int(hit[4]):int(hit[5])]
            hit_header = hit[0]+":"+hit[2]+":"+hit[3]+":"+hit[4]+":"+hit[5]
            gen_seq = seq_dict[gene[1]][int(gene[2]):int(gene[3])]
            gen_header = gene[0]+":"+gene[1]+":"+gene[2]+":"+gene[3]
            outfile.write(">"+str(hit_header)+"\n"+str(hit_seq)+"\n"+">"+str(gen_header)+"\n"+str(gen_seq))
            print(hit_header)
            print(hit_seq)
            print(gen_header)
            print(gen_seq)
