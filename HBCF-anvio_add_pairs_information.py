#!/usr/bin/env python
#This script is designed to add layers to the file that that you created         #
# using the concatenated list of information for genes of the A.vaga genome      #
# created using the HBCF-anvio_assign_variant_positions.py                       #
# This table serves as argument one of the script and argument two is the        #
# name of the file that you want to use for the output.                          #
# so currently the script is run like this                                       #
# python ~/scripts/HBCF-anvio_add_pairs_information.py Layers.txt Layers.add.txt #

# #
# #
# #

import sys

infile = open(sys.argv[1], 'rU') # Name of the table to add information to 
outfile = open(sys.argv[2], 'w') # Name of the outfile that will contain combined information
pairs_info = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles', 'rU') # Table containg additional information
genes_info = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', 'rU')

genes_info_dict = {}
for line in genes_info:
    a = line.strip().split("\t")
    genes_info_dict[a[1]] = a[0:len(a)]

pairs_dict_1 = {} # create a dictionary of the information that you would like to add to the existing table
for line in pairs_info:
    y = line.strip().split("\t")
    pairs_dict_1[y[2]] = y[0:len(y)]

pairs_dict_2 = {}
for line in pairs_info:
    z = line.strip().split("\t")
    pairs_dict_2[z[2]] = z[0:len(z)]

additional_dict = {} # create a dictionary to hold the information that you want to have in the combined table
for line in infile: #search the dictionary containg the additional information and add it the the row
    x = line.strip().split("\t")

    if x[15] == 'intron':
        additional_dict[x[0]] = [x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11], x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],0,0,0,0,0,0,0,0,0,0,0]

    if x[15] not in pairs_dict_1.keys() and x[15] not in pairs_dict_2.keys():
        additional_dict[x[0]] = [x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11], x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],0,0,0,0,0,0,0,0,0,0,0]

    for key in pairs_dict_1:
        if x[15] == key:
            pair_scaf_id = genes_info_dict[pairs_dict_1[key][3]][0]
            additional_dict[x[0]] = [x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11], x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],pairs_dict_1[key][0],pairs_dict_1[key][1],pairs_dict_1[key][2],pairs_dict_1[key][3],pairs_dict_1[key][4], pairs_dict_1[key][5], pairs_dict_1[key][6], pairs_dict_1[key][7], pairs_dict_1[key][8],pairs_dict_1[key][9], pair_scaf_id]
    
    for key in pairs_dict_2:
        if x[15] == key:
            pair_scaf_id = genes_info_dict[pairs_dict_2[key][3]][0]
            additional_dict[x[0]] = [x[0], x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11], x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19],pairs_dict_2[key][0],pairs_dict_2[key][1],pairs_dict_2[key][2],pairs_dict_2[key][3],pairs_dict_2[key][4], pairs_dict_2[key][5], pairs_dict_2[key][6], pairs_dict_2[key][7], pairs_dict_2[key][8],pairs_dict_2[key][9], pair_scaf_id]
    

outfile.write("contig"+"\t"+"split_name"+"\t"+"AV11"+"\t"+"AV12"+"\t"+"AV13"+"\t"+"AV31"+"\t"+"AV32"+"\t"+"AV33"+"\t"+"AV34"+"\t"+"AV35"+"\t"+"Av36"+"\t"+"AVC1"+"\t"+"AVC3"+"\t"+"AVT0"+"\t"+"REF"+"\t"+"CODING"+"\t"+"GENE"+"\t"+"START"+"\t"+"SNP_POS"+"\t"+"STOP"+"\t"+"SNP_CONSENSUS"+"\t"+"block_number"+"\t"+"gene_number"+"\t"+"gene_ID_1"+"\t"+"gene_ID_2"+"\t"+"E-val"+"\t"+"Ka"+"\t"+"Ks"+"\t"+"ratio_Ka-Ks"+"\t"+"ID_nucleotide"+"\t"+"ID_amino"+"\t"+"scaffold_pair"+"\n")#write the header infromation

for key in additional_dict.keys():# write the dictionary to a file
    outfile.write(key+"\t"+str(additional_dict[key][0])+"\t"+str(additional_dict[key][1])+"\t"+str(additional_dict[key][2])+"\t"+str(additional_dict[key][3])+"\t"+str(additional_dict[key][4])+"\t"+str(additional_dict[key][5])+"\t"+str(additional_dict[key][6])+"\t"+str(additional_dict[key][7])+"\t"+str(additional_dict[key][8])+"\t"+str(additional_dict[key][9])+"\t"+str(additional_dict[key][10])+"\t"+str(additional_dict[key][11])+"\t"+str(additional_dict[key][12])+"\t"+str(additional_dict[key][13])+"\t"+str(additional_dict[key][14])+"\t"+str(additional_dict[key][15])+"\t"+str(additional_dict[key][16])+"\t"+str(additional_dict[key][17])+"\t"+str(additional_dict[key][18])+"\t"+str(additional_dict[key][19])+"\t"+str(additional_dict[key][20])+"\t"+str(additional_dict[key][21])+"\t"+str(additional_dict[key][22])+"\t"+str(additional_dict[key][23])+"\t"+str(additional_dict[key][24])+"\t"+str(additional_dict[key][25])+"\t"+str(additional_dict[key][26])+"\t"+str(additional_dict[key][27])+"\t"+str(additional_dict[key][28])+"\t"+str(additional_dict[key][29])+"\t"+str(additional_dict[key][30])+"\n")
