#!/usr/bin/env python

import sys

gene = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', 'rU') #open the genes file
alleles = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles', 'rU') #open the alleles file
ohno = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.ohnologues', 'rU') #open the ohnologues file
unpaired_blast = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/BLAST_RESULTS_UNPAIRED_SCAFFOLDS_1-100.txt', 'rU')

print("creating a beautiful table for you - called PAIRS_CORRECTED.txt")
allele_list = [] # make a list of the alleles
for line in alleles:
    a = line.strip().split("\t")
    allele_list.append([a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]])

ohno_list = [] # make a list of the ohnologues
for line in ohno:
    b = line.strip().split("\t")
    ohno_list.append([b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9]])

gene_dict = {} # create a dictionary for the file containing all of the genes
for line in gene:
    x = line.strip().split("\t")
    #if x[0] == scaffold_genes:
    gene_dict[x[1]] = x[0:len(x)]


unpaired = [] # create a list of unpaired from the BLAST output of unpaired genes (created by HBCF-gene_and_pair_mining.py)
for line in unpaired_blast:
    x = line.strip().split(",")
    header = x[0].split("_") # the header contains the gene name, scaffold ID, and length of gene
    unpaired.append([header[0],header[1]+'_'+header[2],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]])

for line in unpaired:
    print line

List_from_genes = [] # contains the scaffold1, scaffold1, start1, stop1, start2, stop2, e-value, and %Identity


for line in allele_list: # 
    if line[0] in gene_dict.keys():
        List_from_genes.append([gene_dict[line[0]][0],gene_dict[line[1]][0], gene_dict[line[0]][1], gene_dict[line[1]][1], gene_dict[line[0]][2], gene_dict[line[0]][3], gene_dict[line[1]][2], gene_dict[line[1]][3], line[2],line[6]])


for line in List_from_genes:
    if line[0] in gene_dict.keys():
        List_from_genes.append([gene_dict[line[0]][0],gene_dict[line[1]][0], gene_dict[line[0]][1], gene_dict[line[1]][1], gene_dict[line[0]][2], gene_dict[line[0]][3], gene_dict[line[1]][2], gene_dict[line[1]][3], line[2], line[6]])

### I can insert something here to check if the gene overlaps with an existing gene

num = 0
for line in unpaired:
    if line[0] in gene_dict.keys() and int(line[4]) > int(200):
        List_from_genes.append([line[1], line[2], line[0], 'NOPAIR'+'_'+str(num), gene_dict[line[0]][2], gene_dict[line[0]][3], line[9], line[10], line[11], line[3]])
        num += 1

outfile = open('PAIRS_CORRECTED.txt', 'w')
outfile.write("scaffold_1"+"\t"+"scaffold_2"+"\t"+"GENE_1"+"\t"+"GENE_2"+"\t"+"start_1"+"\t"+"stop_1"+"\t"+"start_2"+"\t"+"stop_2"+"\t"+"e-val"+"\t"+"percent_ID"+"\n")
for line in List_from_genes:
    x = int(line[5]) 
    y = int(line[7])
    if x != y:
        outfile.write(str(line[0])+"\t"+str(line[1])+"\t"+str(line[2])+"\t"+str(line[3])+"\t"+str(line[4])+"\t"+str(line[5])+"\t"+str(line[6])+"\t"+str(line[7])+"\t"+str(line[8])+"\t"+str(line[9])+"\n")
