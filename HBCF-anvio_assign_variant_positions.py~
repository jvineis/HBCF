#!/user/bin/env python

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# This script is designed to add an additioal layer to the basic output from the anvio command ##
# anvi-export-table CONTIGS.db --table splits_basic_info -o basic.txt                          ##
# The output table has already been altered to include whether or not the split falls within   ##
# a gene or an intron and the name of the gene using the following scripts.                    ##
# python ~/scripts/au-add_to_basic_info_table.py basic.txt gene_position_list.txt              ##
# We will use the output from this script "genes_per_split.txt" as the input and the matchinig ##
# pairs file that contains the SNPs and their locations within the gene to give a rough        ##
# of where the SNPs fall within each of the splits (if any).  Matching pairs file:             ## 
# /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/TEST_DIR/MATCHING_PAIRS_SNPS.txt   ##
# The script works on an individual gene that you provide and a list of genes can be processed ##
# in a for loop.  I run the command like this, where I have created a list of genes to run     ##
# where each gene is on a single line of a file.                                               ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# for i in `cat GENES.txt`; do python ~/scripts/HBCF-anvio_assign_variant_positions.py genes_per_split.txt $i; done
# Then run just for the introns
# python ~/scripts/HBCF-anvio_assign_variant_positions.py genes_per_split.txt intron
#
# Finally bring them all together 
#cat *snps_per_split.txt intron_snps_per_split.txt > Additional_anvio_layers.txt
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


import sys 

anvio_table = open(sys.argv[1], 'rU')
snp_positions = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/TEST_DIR/MATCHING_PAIRS_SNPS.txt', 'rU')
lookup_gene = sys.argv[2]
anvio_dict = {}

with anvio_table as x:
    for line in x:
        line = line.strip().split("\t")
        anvio_dict[line[0]] = line[1:len(line)] # now we have a dictionary that looks like this
## key [ values ] 
## gnl_Avscaffoldsv2_scaffold_1_1087316bp_split_02897 ['INTRON', 'intron', '126390', '127324'] 

filename = lookup_gene+"_snps_per_split.txt" # create a file name based on the input gene
outfile = open(filename, 'w')
if lookup_gene == "intron":
    n = 0
    for key in anvio_dict.keys():
        if anvio_dict[key][0] == 'INTRON':
        #print (n,key)
            start = int(anvio_dict[key][2])+n
            end = int(anvio_dict[key][3])+n
        #print (start, end)
            outfile.write(str(key) +"\t"+ str(anvio_dict[key][0]) +"\t"+ str(anvio_dict[key][1]) +"\t"+ str(anvio_dict[key][2]) +"\t"+ str(anvio_dict[key][3]) +"\t"+ str(start) +"\t"+"consensus"+"\n")
            n += 1
    genes = [] # a list to hold the genes
    for key in anvio_dict.keys(): #start by going through each of the keys
        if anvio_dict[key][0] == 'GENE': #If the value is a gene, this is great
            x = anvio_dict[key][1] 
            if x in genes:
                next
            else:
                genes.append(x)
else:
    gene_dict = {}
    hit_list = {}
    no_snp_list = {}
    a = 0
    for key in anvio_dict.keys():
        if anvio_dict[key][1] == lookup_gene:
            split_length = int(anvio_dict[key][3])-int(anvio_dict[key][2])
            gene_dict[key] = (anvio_dict[key][1],a, a + split_length)
            a += split_length
                 
    for line in snp_positions: # This is the matching pairs table - it is quite large ~ 450000 lines
        x = line.strip().split("\t") # make the lines parseable
        for key in gene_dict.keys(): # These are the keys in the dictionary - key = split name : values = the GENE name and start stop positions of the gene
            if gene_dict[key][0] == x[0] and int(x[1]) > int(gene_dict[key][1]) and int(x[1]) < int(gene_dict[key][2]): # looking for snps that have the gene name and fall within the start and stop positions of the gene
                hit_list[key] = [str(key), str("EXON"), str(gene_dict[key][0]), str(gene_dict[key][1]), str(x[1]), str(gene_dict[key][2]), str("snp")] # append the goods to the hit list
                print(gene_dict[key], "snp") # print for fun
            #else:
            #    no_snp_list[key] = [str(key), str("EXON"), str(gene_dict[key][0]), str(gene_dict[key][1]), str(x[1]), str(gene_dict[key][2]), str("consensus")]
            #    print (gene_dict[key], "consensus")
    for key in hit_list.keys(): ## write each of the snps to the table with all the good information for anvio :)
        outfile.write(str(key) +"\t"+str(hit_list[key][1]) +"\t"+ str(hit_list[key][2])+"\t"+str(hit_list[key][3]) +"\t"+str(hit_list[key][4]) +"\t"+ str(hit_list[key][5])+"\t"+str(hit_list[key][6]) + "\n")
    for key in gene_dict.keys():
        if key in hit_list.keys():
            next
        else:
            no_snp_list[key] = [str(key), str("EXON"), str(gene_dict[key][0]), str(gene_dict[key][1]), str(x[1]), str(gene_dict[key][2]), str("consensus")]
            print(gene_dict[key], "consensus")
    for key in no_snp_list.keys():
        outfile.write(str(key) +"\t"+str(no_snp_list[key][1]) +"\t"+ str(no_snp_list[key][2])+"\t"+str(no_snp_list[key][3]) +"\t"+str(no_snp_list[key][4]) +"\t"+ str(no_snp_list[key][5])+"\t"+str(no_snp_list[key][6]) + "\n")
