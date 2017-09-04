#!/bin/usr/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = '''Create a bash script to run genewise on a previously unannotated and recently allele paired region of the genome''')
parser.add_argument('fasta_file_list', help='The ALLELES_NEW_GENES.txt or ALLELES_EXISTING_GENES.txt')
parser.add_argument('-dna', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/FASTA/',help ='the full path to the directory containing the fasta files in the fasta_file_list')
parser.add_argument('-pep', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/PEP/', help='the full path to the directory containing the peptide sequence of the known allele')
parser.add_argument('-out', default = '.', help = 'the full path of the directory where you want the genewise shell script to be written')
args = parser.parse_args()

file_names = open(args.fasta_file_list, 'rU')
path_to_seq = args.dna

for name in file_names:
    x = name.strip().split("\t")
    peptide_fasta_name = x[0]+".fa"
    if int(x[4]) < int(x[5]):
        dna_fasta_name = x[2]+"_"+x[3]+"_"+x[4]+"_"+x[5]+".fa"
    if int(x[4]) > int(x[5]):
        dna_fasta_name = x[2]+"_"+x[3]+"_"+x[5]+"_"+x[4]+".fa"
    outfile = open(args.out+"/"+x[2]+".genewise.sh", 'w')
 
    command_line = "#!/bin/bash" +"\n" +"genewise %s/%s %s/%s -gene worm.gf -cdna -pep -both -pretty" %(args.pep, peptide_fasta_name, args.dna, dna_fasta_name) + " > %s.genewise.out.txt" %(dna_fasta_name) 
    outfile.write(str(command_line))
