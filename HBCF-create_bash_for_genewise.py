#!/bin/usr/env python

import sys
from Bio import SeqIO
## Takes as input a table with a single column containing the list fasta files names that 
## you want to run throught genewise.  

## The output is a bash script that you can run throught genewise.

file_names = open(sys.argv[1], 'rU')
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
    
    command_line = "#!/bin/bash" +"\n" +"genewise /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/PEP/%s /workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/FASTA/%s -gene worm.gf -u 1 -v %s -embl -cdna -pep" %(peptide_fasta, name, length[0]) + " > %s.genewise.out.txt" %(name)
    outfile.write(str(command_line))
