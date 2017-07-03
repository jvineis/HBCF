#!/usr/bin/env python
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# I created this script to make a file that can be compatible with anvi-interactive -A.                  ##
# In the current form, it takes two arguments 1. the basic file created from anvio.                      ## 
# anvi-export-table CONTIGS.db --table splits_basic_info -o basic.txt                                    ##
# The file is converted to a dictionary with the individual splits serving as the keys.  I also          ##
# created a file containing the position of genes in my fasta file that I want to add as a layer         ##
# to the display.  I made sure that the coordinates of the genes match the split and not the genome      ##
# details of this script etc can be found here                                                           ##
# /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GC_MAPPING_VISUALIZATION/AVSCAFFOLD_1/README ##
# and here                                                                                               ##
# python ~/scripts/HBCF-genome-filter.py [scaffold_of_interest]                                          ##
#
# If the position of start and stop position of the split falls within the range of the gene coordinates ##
# The gene split name are appended to the output "genes_per_split.txt".  If the split is outside the     ##
# boundaries of all genes, "intron" is appended to the table.

import sys

basic_infile = open(sys.argv[1], 'rU')
geneposition_infile = open(sys.argv[2], 'rU')
outfile = open(sys.argv[3], 'w')

basic_dict = {}
with basic_infile as x:
    firstline = x.readline()
    for line in x:
        x = line.strip().split('\t')
        basic_dict[x[0]] = x[1:int(len(x))]


hits = [] # a list to store the split name and the gene name associated with the split
hitsd = {} # a dictionary to hold the key of the basic_dict and the gene name 
for line in geneposition_infile: # run through the gene positions file containing the gene name, start and stop positoin
    x = line.strip().split("\t") # fix the line into a parsable object
    for key in basic_dict.keys(): # evaluate the values for each key to test whether the split is within the range of the gene
        if int(x[1]) < int(basic_dict[key][1]) and int(x[2]) > int(basic_dict[key][2]): # if the start of the gene x[1] is less than the start of the split basic_dict[key][1]
                                                                                        # and the end of the gene x[2] is more than the end of the split basic_dict[key][2]
            hits.append([key, x[0]]) #append the key and name of the gene to the hit list
            hitsd[key] = x[0] #do the same for a dictionary

for line in hits:# run through the lines of the hits
    x = line# identify the line
    for key in basic_dict.keys(): # if the key is equal to the hit write the information to the genes_per_split table
        if x[0] == key:
            outfile.write(str(x[0])+"\t"+"GENE"+"\t"+str(x[1])+"\t"+str(basic_dict[key][1])+"\t"+str(basic_dict[key][2])+"\n")


for key in basic_dict.keys(): # very simply, if the dictionary key is not associated with a gene. write the details of the split and identify as an intron.
    if key in hitsd.keys():
        print( "found this split", key)
    else:
        print("not fount", key)
        outfile.write(str(key)+"\t"+"INTRON"+"\t"+"intron"+"\t"+str(basic_dict[key][1])+"\t"+str(basic_dict[key][2])+"\n")
