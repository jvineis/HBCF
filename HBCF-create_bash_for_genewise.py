#!/bin/usr/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description = '''Create a bash script to run genewise on a previously unannotated and recently allele paired region of the genome''')
parser.add_argument('fasta_file_list', help='a single column file of fasta file names')
parser.add_argument('path_to_sequences', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/FASTA/',help ='the full path to the directory containing the fasta files in the fasta_file_list')
parser.add_argument('path_to_peptides', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/PEP/', help='the full path to the directory containing the peptide sequence of the known allele')
args = parser.parse_args()

file_names = open(args.fasta_file_list, 'rU')
path_to_seq = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/FASTA/'

for name in file_names:
    x = name.strip().split("_")
    name = name.strip()
    peptide_fasta = x[0][0:14] + "001.fa"
    print(name)
    print (peptide_fasta)
    outfile = open(x[0]+".genewise.sh", 'w')
    length = []
    for rec in SeqIO.parse(open(path_to_seq+name, 'rU'), 'fasta'):
        length.append(len(rec.seq))
 
    command_line = "#!/bin/bash" +"\n" +"genewise %s/%s %s/%s -gene worm.gf -u 1 -v %s -embl -cdna -pep" %(args.path_to_peptides, peptide_fasta, args.path_to_sequences, name, length[0]) + " > %s.genewise.out.txt" %(name) 
#    command_line = "#!/bin/bash" +"\n" +"genewise /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/PEP/%s /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/FASTA/%s -gene worm.gf -u 1 -v %s -embl -cdna -pep" %(peptide_fasta, name, length[0]) + " > %s.genewise.out.txt" %(name)
    outfile.write(str(command_line))
