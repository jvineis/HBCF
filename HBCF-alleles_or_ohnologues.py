#!/bin/usr/env python

import sys

infile = open(sys.argv[1], 'rU') ## This is the ALLELES_NEW_GENES.txt, ALLELES_EXISTING_GENES.txt, OHNO_NEW_GENES.txt, OHNO_EXISTING_GENES.txt etc.. 
single_outfile = open(sys.argv[2], 'w') ## This is the file name that you want to use to export only the genes that have a single allele pair.  
multiple_outfile = open(sys.argv[3], 'w') ## This is the file name that you want to use to expopt only the gnese that have multile alleles -- they are terrible!!!!

n = 0
gene_dict = {}
gene_hits = []
for allele in infile:
    x = allele.strip().split("\t")
    gene_hits.append(x[0])
    gene_dict[n] = x[0:len(x)]
    n += 1



for key in gene_dict.keys():
    if int(gene_hits.count(gene_dict[key][0])) == 1:
        print gene_hits.count(gene_dict[key][0]),gene_dict[key]
        single_outfile.write(str(gene_dict[key][0])+"\t"+str(gene_dict[key][1])+"\t"+str(gene_dict[key][2])+"\t"+str(gene_dict[key][3])+"\t"+str(gene_dict[key][4])+"\t"+str(gene_dict[key][5])+"\t"+str(gene_dict[key][6])+"\t"+str(gene_dict[key][7])+"\t"+str(gene_dict[key][8])+"\n")
        
    else:
        multiple_outfile.write(str(gene_dict[key][0])+"\t"+str(gene_dict[key][1])+"\t"+str(gene_dict[key][2])+"\t"+str(gene_dict[key][3])+"\t"+str(gene_dict[key][4])+"\t"+str(gene_dict[key][5])+"\t"+str(gene_dict[key][6])+"\t"+str(gene_dict[key][7])+"\t"+str(gene_dict[key][8])+"\n")
                    


