#!/usr/bin/env python
# This script is designed to select a designated scaffold from the A.vaga genome and write a file that records
# the gene and start stop postion from the "genes.all" file located here
# /groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFile
#  run like this
# python ~/scripts/HBCF-genome-filter.py [scaffold_name e.g. av1]


import sys

#infile = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', 'rU')
infile = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_GENES_AND_PAIRS/GENES.ALL_JHV.txt')

scaffold = sys.argv[1]
outname = scaffold+'_gene_position_list.txt'
indict = {}
for line in infile:
    x = line.strip().split('\t')
    indict[x[0]] = x[0:int(len(x))]
outfile = open(outname, 'w')
for key in indict.keys():
    if indict[key][1] == scaffold: #'av1' and int(indict[key][2]) > 52845 and int(indict[key][3]) < 180171:
        start = int(indict[key][2]) #- 52845
        end = int(indict[key][3]) #- 52845
        outfile.write(str(key) + "\t" + str(start) + "\t" +  str(end) + "\n")

#for key in indict.keys():
#    if indict[key][0] == 'av1' and int(indict[key][2]) > 52845 and int(indict[key][3]) < 180171:
#        start = int(indict[key][2]- 52845)
#        end = int(indict[key][3]- 52845)
#        outfile.write(str(key) + "\t" + str(start) + "\t" +  str(end) + "\n")
