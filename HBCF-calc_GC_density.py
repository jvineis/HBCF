#!/usr/bin/env python

# The purpose of this script is to go through the gene annotation that was created from HBCF-anvio_assign_variant_positions.py ##
# and concatenated into a table of gene and intron positions 'ALL_intron_and_snps_per_split.txt) to determine how many         ## 
# consecutive genes have zero snps                                                                                             ##
# 20 is the maximum number that this scritp can handle currently.                                                              ## 
# The script is run like this                                                                                                  ##
# python ~/scripts/HBCF-calc_GC_density.py GENES.txt ALL_intron_and_snps_per_split.txt                                         ##


import sys

genes = open(sys.argv[1], 'rU') # a list of the genes that are in the scaffold 
snps_per_split = open(sys.argv[2], 'rU') # The splits file contiaining information for each split, whether it contains a SNP or not.  This is a rough
outfile_prefix = sys.argv[3]
#estimate for the number of SNPs because there could be more than one snp in the 25bp split.

gene_list = [] # convert the genes into a list
for gene in genes:
    x = gene.strip().split("\t")
    gene_list.append(x)


gene_dict = {}
count = 0
gene_list = sorted(gene_list)


for gene in gene_list: # Here is how we make a dictionary for each gene with the ten genes downstream as the values 
    print (count, gene)
    if count >= 0 and count <= len(gene_list)-20:
        gene_dict[gene[0]] = gene_list[count][0], gene_list[count+1][0], gene_list[count+2][0], gene_list[count+3][0], gene_list[count+4][0], gene_list[count+5][0], gene_list[count+6][0], gene_list[count+6][0], gene_list[count+7][0], gene_list[count+8][0], gene_list[count+9][0], gene_list[count+10][0], gene_list[count+11][0], gene_list[count+12][0], gene_list[count+13][0], gene_list[count+14][0], gene_list[count+15][0], gene_list[count+16][0], gene_list[count+17][0], gene_list[count+18][0], gene_list[count+19][0]#, gene_list[count+20][0]
    count += 1
    if count > len(gene_list)-20:
        print("this one", gene)

for gene in gene_dict.keys():
    print( gene, gene_dict[gene])

splits_dict = {} # convert the splits infromation into a dictionary, split is the key and all other information is retained as values
for split in snps_per_split:
    x = split.strip().split("\t")
    splits_dict[x[0]] = x[0:len(x)]


snps = [] # make a list for all the gene names whose split is identified as containing a snp.  This is a redundant list indicative of the number of snps in a gene
for key in splits_dict.keys():
    if splits_dict[key][6] == 'snp':
        snps.append(splits_dict[key][2])

snp_counts = {} # Now go through each of the genes in our gene list for the scaffold and find out how many SNPs are identified (at least within 25bp splits) in each of the 10 genes downstream. 
for gene in gene_dict.keys():
    snp_counts[gene] = [gene,snps.count(gene_dict[gene][0]), snps.count(gene_dict[gene][1]), snps.count(gene_dict[gene][2]), snps.count(gene_dict[gene][3]), snps.count(gene_dict[gene][4]), snps.count(gene_dict[gene][5]), snps.count(gene_dict[gene][6]), snps.count(gene_dict[gene][7]), snps.count(gene_dict[gene][8]), snps.count(gene_dict[gene][9]), snps.count(gene_dict[gene][10]), snps.count(gene_dict[gene][11]), snps.count(gene_dict[gene][12]), snps.count(gene_dict[gene][13]), snps.count(gene_dict[gene][14]), snps.count(gene_dict[gene][15]), snps.count(gene_dict[gene][16]), snps.count(gene_dict[gene][17]), snps.count(gene_dict[gene][18]), snps.count(gene_dict[gene][19]), snps.count(gene_dict[gene][20])]

count1 = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
count6 = 0
count7 = 0
count8 = 0
count9 = 0
count10 = 0
count11 = 0
count12 = 0
count13 = 0
count14 = 0
count15 = 0
count16 = 0
count17 = 0
count18 = 0 
count19 = 0
count20 = 0
for snp in snp_counts.keys():
    print (snp_counts[snp])
    if snp_counts[snp][1] == 0 and snp_counts[snp][2] > 0:
        print("a single flanker", snp, snp_counts[snp])
        count1 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] > 0:
        print("a double flanker", snp, snp_counts[snp])
        count2 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] > 0:
        print("a triple flanker", snp, snp_counts[snp])
        count3 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] > 0:
        print("a four flanker", snp, snp_counts[snp])
        count4 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] > 0:
        print("a 5 flanker",snp, snp_counts[snp])
        count5 += 1 
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] > 0:
        print("a 6 flanker",snp, snp_counts[snp])
        count6 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] > 0:
        print("a 7 flanker",snp, snp_counts[snp])
        count7 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] > 0:
        print("a 8 flanker",snp, snp_counts[snp])
        count8 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] > 0:
        print("a 9 flanker",snp, snp_counts[snp])
        count9 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] > 0:
        print("a 10 flanker",snp, snp_counts[snp])
        count10 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] > 0:
        print("a 11 flanker",snp, snp_counts[snp])
        count11 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] > 0:
        print("a 12 flanker",snp, snp_counts[snp])
        count12 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] ==  0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] > 0:
        print("a 13 flanker",snp, snp_counts[snp])
        count13 += 1    
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] ==0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] > 0:
        print("a 14 flanker",snp, snp_counts[snp])
        count14 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] ==0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] > 0:
        print("a 15 flanker",snp, snp_counts[snp])
        count15 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] == 0 and snp_counts[snp][17] > 0:
        print("a 16 flanker",snp, snp_counts[snp])
        count16 += 1   
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] == 0 and snp_counts[snp][17] == 0 and snp_counts[snp][18] > 0:
        print("a 17 flanker",snp, snp_counts[snp])    
        count17 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] == 0 and snp_counts[snp][17] == 0 and snp_counts[snp][18] == 0 and snp_counts[snp][19] > 0:
        print("a 18 flanker",snp, snp_counts[snp])
        count18 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] == 0 and snp_counts[snp][17] == 0 and snp_counts[snp][18] == 0 and snp_counts[snp][19] == 0 and snp_counts[snp][20] > 0:
        print("a 19 flanker",snp, snp_counts[snp])
        count19 += 1
    elif snp_counts[snp][1] == 0 and snp_counts[snp][2] == 0 and snp_counts[snp][3] == 0 and snp_counts[snp][4] == 0 and snp_counts[snp][5] == 0 and snp_counts[snp][6] == 0 and snp_counts[snp][7] == 0 and snp_counts[snp][8] == 0 and snp_counts[snp][9] == 0 and snp_counts[snp][10] == 0 and snp_counts[snp][11] == 0 and snp_counts[snp][12] == 0 and snp_counts[snp][13] == 0 and snp_counts[snp][14] == 0 and snp_counts[snp][15] == 0 and snp_counts[snp][16] == 0 and snp_counts[snp][17] == 0 and snp_counts[snp][18] == 0 and snp_counts[snp][19] == 0 and snp_counts[snp][20] == 0 and snp_counts[snp][21] > 0:
        print("a 20 flanker",snp, snp_counts[snp])
        count20 += 1

print(count20,count19,count18,count17,count16,count15,count14,count13,count12,count11,count10,count9,count8,count7,count6,count5,count4,count3,count2,count1)

real20 = count20
if count19 >= real20 and count19 != 0:
    real19 = count19 - count20
#elif count19 <= 1 and real20 > count19:
#    real19 = 0
else:
    real19 = count19

if count18 >= real19 and count18 != 0:
    real18 = count18 - count19
#elif count18 <= 1 and real19 > count18:
#    real18 = 0
else:
    real18 = count18

if count17 >= real18 and count17 != 0:
    real17 = count17 - count18
#elif count17 <= 1 and real18 > count17:
#    real17 = 0
else:
    real17 = count17

if count16 >= real17 and count16 != 0:
    real16 = count16 - count17
#elif count16 <= 1 and real17 > count16:
#    real16 = 0
else:
    real16 = count16

if count15 >= real16 and count15 != 0:
    real15 = count15 - count16
#elif count15 <= 1 and real16 > count15:
#    real15 = 0
else:
    real15 = count15

if count14 >= real15 and count14 != 0:
    real14 = count14 - count15
#elif count14 <= 1 and real15 > count14:
#    real14 = 0
else:
    real14 = count14

if count13 >= real14 and count13 != 0:
    real13 = count13 - count14
#elif count13 <= 1 and real14 > count13:
#    real13 = 0
else:
    real13 = count13

if count12 >= real13 and count12 != 0:
    real12 = count12 - count13
#elif count12 <= 1 and real13 > count12:
#    real12 = 0
else:
    real12 = count12

if count11 >= real12 and count11 != 0:
    real11 = count11 - count12
#elif count11 <= 1 and real12 > count11:
#    real11 = 0
else:
    real11 = count11

if count10 >= real11 and count10 != 0:
    real10 = count10 - count11
#elif count10 <= 1 and real11 > count10:
#    real10 = 0
else:
    real10 = count10

if count9 >= real10 and count9 != 0:
    real9 = count9 - count10
#elif count9 <= 1 and real10 > count9:
#    real9 = 0
else:
    real9 = count9

if count8 >= real9 and count8 != 0:
    real8 = count8 - count9
#elif count8 <= 1 and real9 > count8:
#    real8 = 0
else:
    real8 = count8

if count7 >= real8 and count7 != 0:
#    print("count7 is %d and real8 is %d" %(count7, real8)) 
    real7 = count7 - count8
#elif count7 <= 1 and real8 > count7:
#    real7 = 0
else:
    real7 = count7

if count6 >= real7 and count6 != 0:
    real6 = count6 - count7
#elif count6 <= 1 and real7 > count6:
#    real6 = 0
else:
    real6 = count6

if count5 >= real6 and count5 != 0:
    real5 = count5 - count6
#elif count5 <= 1 and real6 > count5:
#    real5 = 0
else:
    real5 = count5

if count4 >= real5 and count4 != 0:
    real4 = count4 - count5
#elif count4 <= 1 and real5 > count4:
#    real4 = 0
else:
    real4 = count4

if count3 >= real4 and count3 != 0:
    real3 = count3 - count4
#elif count3 <= 1 and real4 > count3:
#    real3 = 0
else:
    real3 = count3

if count2 >= real3 and count2 != 0:
    real2 = count2 - count3
#elif count2 <= 1 and real3 > count2:
#    real2 = 0
else:
    real2 = count2

if count1 >= real2 and count1 != 0:
    real1 = count1 - count2
#elif count1 <= 1 and real2 > count1:
#    real1 = 0
else:
    real1 = count1




print(real20, real19, real18, real17, real16, real15, real14, real13, real12,real11, real10,real9, real8, real7, real6, real5, real4, real3, real2, real1)

outfile = open(outfile_prefix+'_GC_COUNTS.txt', 'w')
outfile.write('20-GENE'+'\t'+'19-GENE'+'\t'+'18-GENE'+'\t'+'17-GENE'+'\t'+'16-GENE'+'\t'+'15-GENE'+'\t'+'14-GENE'+'\t'+'13-GENE'+'\t'+'12-GENE'+'11-GENE'+'\t'+'10-GENE'+'\t'+'9-GENE'+'\t'+'8-GENE'+'\t'+'7-GENE'+'\t'+'6-GENE'+'\t'+'5-GENE'+'\t'+'4-GENE'+'\t'+'3-GENE'+'\t'+'2-GENE'+'\t'+'1-GENE'+'\n')
outfile.write(outfile_prefix+'\t'+str(real20)+'\t'+str(real19)+'\t'+str(real18)+'\t'+str(real17)+'\t'+str(real16)+'\t'+str(real15)+'\t'+str(real14)+'\t'+str(real13)+'\t'+str(real12)+'\t'+str(real11)+'\t'+str(real10)+'\t'+str(real9)+'\t'+str(real8)+'\t'+str(real7)+'\t'+str(real6)+'\t'+str(real5)+'\t'+str(real4)+'\t'+str(real3)+'\t'+str(real2)+'\t'+str(real1)+'\n') 






