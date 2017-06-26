#!/usr/bin/env python

#This script selects genes that are not found in the A.vaga genes or ohnologues file and creates a fasta file#
#for the individual gene to be run through BLAST.  The user selects that target scaffold number and runs it  #
# like this         
# python HBCF-gene_and_pair_mining.py 1 
# will return the results for all genes in scaffold 1 that are the targets of the search.

import sys
from Bio import SeqIO

gene = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', 'rU') #open the genes file
alleles = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles', 'rU') #open the alleles file
ohno = open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.ohnologues', 'rU') #open the ohnologues file
scaffold = open('/groups/rotifer/Avgenome/Genoscope/v2/Joes_Adineta_vaga_v2.0.scaffolds.fa', 'rU') #open the fasta file

scaffold_name = sys.argv[1] # this should just be a single numeric value 1,2,3.. etc.
scaffold_genes = 'av'+scaffold_name
scaffold_fasta = 'scaffold_'+scaffold_name+'_'
print ("searching for lonely pairs in %s" %(scaffold_fasta))

allele_list = [] # make a list of the alleles
for line in alleles:
    a = line.strip().split("\t")
    allele_list.append(a[2])
    allele_list.append(a[3])

ohno_list = [] # make a list of the ohnologues
for line in ohno:
    b = line.strip().split("\t")
    ohno_list.append(b[2])
    ohno_list.append(b[3])

gene_dict = {} # create a dictionary for the file containing all of the genes
for line in gene:
    x = line.strip().split("\t")
    if x[0] == scaffold_genes:
        gene_dict[x[1]] = x[0:len(x)]


## Now we want to search our genes and determine whether the are in the allele list and the ohno_list or unpaired

unpaired = []
for key in gene_dict.keys():
    if gene_dict[key][1] in allele_list:
        next
    elif gene_dict[key][1] in ohno_list:
        next
    else:
        unpaired.append([gene_dict[key][1],gene_dict[key][2],gene_dict[key][3]])

target_seq = [] ## Select the target sequence from our geneome
for record in SeqIO.parse(scaffold, "fasta"):
    if scaffold_fasta in record.id:
        target_seq.append(record.seq)

def select_from_fasta(name, gene, start, end): ## A function to select the region of sequence from the genome that corresponds with the scaffold.
    outfile_name = gene+'_'+name
    outfile = open(outfile_name, 'w')
    length = int(end)-int(start)
    x = target_seq[0][int(start):int(end)]
    outfile.write(">"+gene+'_'+name+str(length)+'bp'+"\n"+str(x))
    return outfile

for line in unpaired: ## Run each gene that is not found in either pairs file throught the function to create and individual fasta file.
    select_from_fasta(scaffold_fasta, line[0], line[1], line[2])

