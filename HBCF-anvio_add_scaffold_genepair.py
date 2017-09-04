#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='''Adding the scaffold and allele pair information for each split in the display''')
parser.add_argument('anvio_vis_table', help = 'the Anvio_display_layers.txt file usually found in the visualization directory for each scaffold')
parser.add_argument('-pairs_table', help = 'the table with gene1, scaffold1, start1, end1, gene2, scaffold2, start2, end2 - usually created by ~/scripts/HBCF-combine_table_info.py', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/gene_pairs_with_scaffold_info.txt')
parser.add_argument('-outfile', help = 'the file to write the table with the best anvio information ever', default = 'Anvio_display_layers_with_scaffolds.txt')
parser.add_argument('-alien', help = 'the table containing the alien index for each gene', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/AI.tab') 
args = parser.parse_args()

pairs_dict = {}
pairs_table = open(args.pairs_table, 'rU')
count = 0
pairs_list = []
for line in pairs_table:
    x = line.strip().split('\t')
    pairs_dict[count] = x[0:len(x)]
    pairs_list.append(x[0])
    pairs_list.append(x[4])
    count += 1

outfile = open(args.outfile, 'w')
outfile.write("split"+'\t'+"gene"+'\t'+"SNP"+'\t'+"REF"+'\t'+"HER"+'\t'+"REF_trans"+'\t'+"pair_a"+'\t'+"scaffold_a"+'\t'+"start_a"+'\t'+"end_a"+'\t'+"pair_b"+'\t'+"scaffold_b"+'\t'+"start_b"+'\t'+"end_b"+'\n')
split_table = open(args.anvio_vis_table, 'rU')


for line in split_table:
    x = line.strip().split('\t')
    if x[1] == "INTRON":
        outfile.write(x[0]+'\t'+x[2]+'\t'+x[6]+'\t'+x[7]+'\t'+x[8]+'\t'+x[9]+'\t'+"intergenic"+'\t'+"unknown"+'\t'+"0"+'\t'+"0"+'\t'+"intergenic"+'\t'+"unknown"+'\t'+"0"+'\t'+"0"+'\n')
        
    if x[2] not in pairs_list and x[1] == "EXON":
        print(x[2])
        outfile.write(x[0]+'\t'+x[2]+'\t'+x[6]+'\t'+x[7]+'\t'+x[8]+'\t'+x[9]+'\t'+"unpaired"+'\t'+"unknown"+'\t'+"0"+'\t'+"0"+'\t'+"unpaired"+'\t'+"unknown"+'\t'+"0"+'\t'+"0"+'\n')

    if x[2] in pairs_list:
        split_list = []
        for key in pairs_dict.keys():
            if x[2] == pairs_dict[key][0] and x[0] not in split_list:
                outfile.write(x[0]+'\t'+x[2]+'\t'+x[6]+'\t'+x[7]+'\t'+x[8]+'\t'+x[9]+'\t'+pairs_dict[key][0]+'\t'+pairs_dict[key][1]+'\t'+str(pairs_dict[key][2])+'\t'+str(pairs_dict[key][3])+'\t'+pairs_dict[key][4]+'\t'+pairs_dict[key][5]+'\t'+str(pairs_dict[key][6])+'\t'+str(pairs_dict[key][7])+'\n')
                split_list.append(x[0])
            if x[2] == pairs_dict[key][4] and x[0] not in split_list:
                outfile.write(x[0]+'\t'+x[2]+'\t'+x[6]+'\t'+x[7]+'\t'+x[8]+'\t'+x[9]+'\t'+pairs_dict[key][4]+'\t'+pairs_dict[key][5]+'\t'+str(pairs_dict[key][6])+'\t'+str(pairs_dict[key][7])+'\t'+pairs_dict[key][0]+'\t'+pairs_dict[key][1]+'\t'+str(pairs_dict[key][2])+'\t'+str(pairs_dict[key][3])+'\n')
                split_list.append(x[0])
            else:
                next



