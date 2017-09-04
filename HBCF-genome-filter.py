#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='''This script selects a designated scaffold from the A.vaga genome and writes a file that records the start and stop position from the "genes.all" file located in a defined file such as the genes.all file.  The outfile produced includes the gene name, start of the gene and end of the gene.  The file is named after the scaffold that you are searching -eg. av1 will yield an output called av1_gene_position_list.txt''')
parser.add_argument('--genes', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/genes.all', help = 'This is a file formatted in a similar way to the genes.all file here - "/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/genes.all". Its a four column txt file containing the scaffold_id, Gene_name, gene_start, gene_stop')
parser.add_argument('scaffold', help = 'The name of the scaffold in the A. vaga genome that you would like to analyze. eg "av1"')
args = parser.parse_args()

infile = open(args.genes, 'rU') # the genes file with start and stop positions
scaffold = args.scaffold # the name of the scaffold of interest
outname = scaffold+'_gene_position_list.txt'
indict = {}

for line in infile: 
    x = line.strip().split('\t')
    print(x)
    indict[x[1]] = x[0:int(len(x))]

outfile = open(outname, 'w')

for key in indict.keys():
    if indict[key][0] == scaffold: #'av1' and int(indict[key][2]) > 52845 and int(indict[key][3]) < 180171:
        start = int(indict[key][2]) #- 52845
        end = int(indict[key][3]) #- 52845
        outfile.write(str(key) + "\t" + str(start) + "\t" +  str(end) + "\n")
