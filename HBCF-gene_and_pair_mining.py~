#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Export fasta files for unpaired genes in a defined scaffold of the A.vaga genome''')
parser.add_argument('scaffold', type=str, help = 'the scaffold number of the A.vaga genome that you want to analyze, e.g. 1')
parser.add_argument('--genes_table', default='/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', help = 'The full path to the genes.all file of the A.vaga genome v2 or a tab separated table: first column [scaffold id]; second column [gene_id]; third column [gene start position]; fourth column [gene end position] e.g. [av1][GSADVT00000034001][123][225]')
parser.add_argument('--alleles', default='/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles', help = 'The table of allele pairs from the v2 Consortium files named pairs.alleles or a table: first column [random]; second column[random]; third column [gene_1]; fourth column [gene_2] e.g. [0][0][GSADVT00000227001][GSADVT00022299001]')
parser.add_argument('--ohno', default='/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.ohnologues', help = 'The table of the ohonologues from the v2 Consortium files named pairs.ohnologues, same file format as the alleles table')
parser.add_argument('--genome', default='/groups/rotifer/Avgenome/Genoscope/v2/Joes_Adineta_vaga_v2.0.scaffolds.fa', help = 'The v2 A.vaga genome fasta file or a fasta file with headers named similar to this >scaffold_1_1087316_bp.')
args = parser.parse_args()

gene = open(args.genes_table, 'rU') #open the genes file
alleles = open(args.alleles, 'rU') #open the alleles file
ohno = open(args.ohno, 'rU') #open the ohnologues file
scaffold = open(args.genome, 'rU') #open the fasta file

scaffold_name = args.scaffold# this should just be a single numeric value 1,2,3.. etc.
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
