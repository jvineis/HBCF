#!/bin/usr/env python

import argparse

parser = argparse.ArgumentParser(description='''Separates the alleles into those that obey the tetraploid rule and those that are hemizygous''')
parser.add_argument('-directory', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/', help = 'The directory containing the tables called "ALLELES_EXISTING_GENES.txt", "ALLELES_NEW_GENES.txt", "OHNO_EXISTING_GENES.txt", "OHNO_NEW_GENES.txt" produced by HBCF-filter_blast_for_pair_correcting.py')
parser.add_argument('-genes', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', help = 'The A.vaga genome consortium file of all identified genes')
parser.add_argument('-output_directory', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/', help = 'The directory where you want your tables to go')
parser.add_argument('-pairs', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.all', help = 'The A.vaga genome consortium file of all pairs')
parser.add_argument('-BLAST', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/BLAST_ALL_SCAFFOLDS_ABOVE_10kb.txt', help = 'the BLAST table containing the hits to unpaired genes, used to double check the results of genes identified a hemizygous')
parser.add_argument('-unpaired_genes', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/UNPAIRED_GENES.txt', help = 'The list of genes that are found on the scaffolds under examination (usually the ones greater than 10kbp) - a single column on data that must look like something like this GSADVT00000301001_scaffold_1_') 
args = parser.parse_args()

new_alleles = open(args.directory + 'ALLELES_NEW_GENES.txt', 'rU')
new_ohno = open(args.directory + 'OHNO_NEW_GENES.txt', 'rU')
ex_alleles = open(args.directory + 'ALLELES_EXISTING_GENES.txt', 'rU')
ex_ohno = open(args.directory + 'OHNO_EXISTING_GENES.txt', 'rU')
pairs = open(args.pairs, 'rU')
genes = open(args.genes, 'rU')
blast = open(args.BLAST, 'rU')
long_scaffold_genes = open(args.unpaired_genes, 'rU')
# A dictionary for the BLAST table to check hemizygous alleles against:
blast_dict = {}
for line in blast: # parse the blast file
    x= line.strip().split(",")
    header = x[0].split("_")
    #blast.append([header[0],header[1]+'_'+header[2], header[3].strip("bp"),  x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8],x[9],x[10],x[11]])
    blast_dict[header[0]+'_'+header[2]+'_'+x[8]+'_'+x[9]] = ['av'+header[2], header[0], x[8], x[9], x[1], x[2], x[3]]

# These are new positions in the genome that passed through the HBCF-filter
# Here we turn all tables into an iterable dictionaries - keys look like this - GSADVT00021595001_96966_97185
alleles_dict = {}
for line in new_alleles:
    na = line.strip('\n').split('\t')
    alleles_dict[na[0]+'_'+na[4]+'_'+na[5]] = na[0:len(na)]
ohno_dict = {}
for line in new_ohno:
    no = line.strip('\n').split('\t')
    ohno_dict[no[0]+'_'+no[4]+'_'+no[5]] =no[0:len(no)]
alleles_exist = {}
for line in ex_alleles:
    na = line.strip('\n').split('\t')
    alleles_exist[na[0]+'_'+na[6]+'_'+na[7]] = na[0:len(na)]
ohno_exist = {}
for line in ex_ohno:
    no = line.strip('\n').split('\t')
    ohno_exist[no[0]+'_'+no[6]+'_'+no[7]] =no[0:len(no)]

print("you have %d new ohnologues" %(len(ohno_dict)))
print("you have %d new alleles" %(len(alleles_dict)))
print("you have %d existing ohnologues" %(len(ohno_exist)))
print("you have %d existing alleles" %(len(alleles_exist)))

# Place all the gene names that found alleles and ohnologues into either an allele or ohnologue name list !! THIS IS A REPETITIVE LIST
found_new_alleles = []
for key in alleles_dict.keys():
    found_new_alleles.append(alleles_dict[key][0])
for key in alleles_exist.keys():
    found_new_alleles.append(alleles_exist[key][0])
print("the sum of all alleles is %d" %(len(found_new_alleles)))
found_new_ohnos = []
for key in ohno_dict.keys():
    found_new_ohnos.append(ohno_dict[key][0])
for key in ohno_exist.keys():
    found_new_ohnos.append(ohno_exist[key][0])
print("the sum of all ohnos is %d"%(len(found_new_ohnos)))


# A dictionary to hold the information in the genes.all table created by th A.vaga Consortium
gene_dict = {}
for line in genes:
    line = line.strip('\n').split('\t')
    k = line[1]+'_'+line[2]+'_'+line[3]
    gene_dict[k] = line[0:len(line)]

genes_in_pairs = [] # A list of the genes in pairs
for line in pairs:
    p = line.strip('\n').split('\t')
    genes_in_pairs.append(p[2])
    genes_in_pairs.append(p[3])

genes_in_large_scaffolds = [] # We only examined the genes that were on scaffolds greater than 10Kbp so we need to filter out genes on all other scaffolds. Otherwise this code will report the genes on smaller scaffolds as hemizygous.  This would be terribly wrong.
for line in long_scaffold_genes:
    ls = line.strip('\n').split('_')
    genes_in_large_scaffolds.append(ls[0])

print ("the total genes in the long scaffolds analyzed is %d"%(len(genes_in_large_scaffolds)))
unpaired = [] # Iterate through the genes and identify the genes that don't have a pair in the A.vaga Consortium "pairs.all" file
for key in gene_dict.keys():
    if gene_dict[key][1] in genes_in_pairs:
        next
    elif gene_dict[key][1] in genes_in_large_scaffolds:
        unpaired.append(gene_dict[key][1])

hemizygous_list = [] # Iterate through the unpaired genes and remove the genes that were identified as ohnonlogues and alleles and show genes that are in the BLAST hit but not identified as alleles or ohnologues
for i in unpaired:
    if i not in found_new_alleles and i not in found_new_ohnos:
        hemizygous_list.append(i)
    
print("found %d hemizygous"%(len(hemizygous_list)))

hemizygous_tracker = []
hemizygotes = open(args.output_directory+'hemizygotes.txt', 'w')
for key in gene_dict.keys():
    if gene_dict[key][1] in hemizygous_list:
        hemizygotes.write(gene_dict[key][0]+'\t'+ gene_dict[key][1]+'\t'+ gene_dict[key][2]+'\t'+ gene_dict[key][3]+'\t'+ gene_dict[key][4]+'\t'+ gene_dict[key][5]+'\t'+ gene_dict[key][6]+'\t'+ gene_dict[key][7]+'\t'+ gene_dict[key][8]+'\t'+ gene_dict[key][9]+'\t'+ gene_dict[key][10]+'\n')
        hemizygous_tracker.append(key)

print("I wrote %d hemizygous genes to your file"%(len(hemizygous_tracker)))

# make a dictionary of all the genes that have a new alleleic hit and record the number of allelic hits
new_allele_counts = {}
n = 0
for key in gene_dict.keys():# This captures all genes that have a new allelic hit and the number of alleleic hits as well as all the information pertaining to the new region identified in the genome
    if gene_dict[key][1] in found_new_alleles:
        new_allele_counts[n] = gene_dict[key][0], gene_dict[key][1], gene_dict[key][2], gene_dict[key][3], gene_dict[key][4], gene_dict[key][5], gene_dict[key][6], gene_dict[key][7], gene_dict[key][8], gene_dict[key][9], gene_dict[key][10], found_new_alleles.count(gene_dict[key][1])
    n += 1

# new_allele_counts looks like this -  8111: ('av194', 'GSADVT00037132001', '299973', '301014', '-19', '53.8426', '27.6', '5', '1', '0', '786', 1)

perfect_tets = open(args.output_directory+'perfect_tets.txt', 'w')
perfect_tets_list = []
for key in new_allele_counts.keys(): # This captures the genes that have one alleleic hit and two ohnologues
    if new_allele_counts[key][1] in found_new_ohnos and new_allele_counts[key][11] == 1:
        a = found_new_ohnos.count(new_allele_counts[key][1])
        if a == 2:
            perfect_tets_list.append(new_allele_counts[key][1])

for key in alleles_exist.keys():
    if alleles_exist[key][0] in perfect_tets_list:
        perfect_tets.write(alleles_exist[key][0]+'\t'+alleles_exist[key][1]+'\t'+alleles_exist[key][2]+'\t'+alleles_exist[key][3]+'\t'+alleles_exist[key][6]+'\t'+alleles_exist[key][7]+'\n')

for key in ohno_exist.keys():
    if ohno_exist[key][0] in perfect_tets_list:
        perfect_tets.write(ohno_exist[key][0]+'\t'+ohno_exist[key][1]+'\t'+ohno_exist[key][2]+'\t'+ohno_exist[key][3]+'\t'+ohno_exist[key][6]+'\t'+ohno_exist[key][7]+'\n')

for key in alleles_dict.keys():
    if alleles_dict[key][0] in perfect_tets_list:
        perfect_tets.write(alleles_dict[key][0]+'\t'+alleles_dict[key][1]+'\t'+alleles_dict[key][2]+'\t'+alleles_dict[key][3]+'\t'+alleles_dict[key][4]+'\t'+alleles_dict[key][5]+'\n')

for key in ohno_dict.keys():
    if ohno_dict[key][0] in perfect_tets_list:
        perfect_tets.write(ohno_dict[key][0]+'\t'+ohno_dict[key][1]+'\t'+ohno_dict[key][2]+'\t'+ohno_dict[key][3]+'\t'+ohno_dict[key][4]+'\t'+ohno_dict[key][5]+'\n')
           

