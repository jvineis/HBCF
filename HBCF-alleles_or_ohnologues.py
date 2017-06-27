#!/bin/usr/env python

import sys
import argparse
parser = argparse.ArgumentParser(description='''Separates the quality blast hits filtered by HBCF-filter_blast_for_pair_correcting.py into potential Ohnologues and Alleles''')
parser.add_argument('allele_table', help = 'one of the four allele or ohnologue tables produced by HBCF-filter_blast_for_pair_correcting.py')
parser.add_argument('single_hit_table', help = 'the table of hits that have a match to only one other postion in the genome')
parser.add_argument('multiple_hit_table', help = 'the table of hits that have a match to multiple position in the genome')
parser.add_argument('--number_of_hits', default = 1, help = 'change this value to 2 when looking for ohnologues')
args = parser.parse_args()

infile = open(args.allele_table, 'rU') ## This is the ALLELES_NEW_GENES.txt, ALLELES_EXISTING_GENES.txt, OHNO_NEW_GENES.txt, OHNO_EXISTING_GENES.txt etc.. 
single_outfile = open(args.single_hit_table, 'w') ## This is the file name that you want to use to export only the genes that have a single allele pair.  
multiple_outfile = open(args.multiple_hit_table, 'w') ## This is the file name that you want to use to expopt only the gnese that have multile alleles -- they are terrible!!!!

n = 0
gene_dict = {}
gene_hits = []
for allele in infile:
    x = allele.strip().split("\t")
    gene_hits.append(x[0])
    gene_dict[n] = x[0:len(x)]
    n += 1



for key in gene_dict.keys():
    if int(gene_hits.count(gene_dict[key][0])) == args.number_of_hits:
        print gene_hits.count(gene_dict[key][0]),gene_dict[key]
        single_outfile.write(str(gene_dict[key][0])+"\t"+str(gene_dict[key][1])+"\t"+str(gene_dict[key][2])+"\t"+str(gene_dict[key][3])+"\t"+str(gene_dict[key][4])+"\t"+str(gene_dict[key][5])+"\t"+str(gene_dict[key][6])+"\t"+str(gene_dict[key][7])+"\t"+str(gene_dict[key][8])+"\n")
        
    else:
        multiple_outfile.write(str(gene_dict[key][0])+"\t"+str(gene_dict[key][1])+"\t"+str(gene_dict[key][2])+"\t"+str(gene_dict[key][3])+"\t"+str(gene_dict[key][4])+"\t"+str(gene_dict[key][5])+"\t"+str(gene_dict[key][6])+"\t"+str(gene_dict[key][7])+"\t"+str(gene_dict[key][8])+"\n")
                    


