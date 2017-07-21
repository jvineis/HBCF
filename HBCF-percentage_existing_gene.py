#!/usr/bin/env python 

import argparse

parser = argparse.ArgumentParser(description = '''Calculate the percentage of ovelap between the genes identified as "new allele pairs" and the existing genes.  The code takes the output of HBCF-filter_blast_for_pair_correcting.py.  Either the the ALLELES_EXISTING_GENES.txt or OHNO_EXISTING_GENES.txt''')
parser.add_argument('ohno_allele_table', help = 'This is this is the greatest freaking table with columns containing [unpaired_gene] [unpaired_scaffold] [new_gene_name] [existing_gene_scaffold] [existing_gene_start][existing_gene_stop] [new_gene_start] [new_gene_stop] [%identity]')
parser.add_argument('out_table', help = 'The name of the table where you want to write the data')

args = parser.parse_args()
table = open(args.ohno_allele_table, 'rU')
output = open(args.out_table, 'w')
output.write("unpaired_gene"+'\t'+"existing_gene"+'\t'+"percent_overalap"+'\n')
with table as f:
    firstline = f.readline()
    for line in table:
        b = line.strip('\n').split('\t')
        if int(b[4]) < int(b[5]):
            bases_a = range(int(b[4]),int(b[5]))
        if int(b[4]) > int(b[5]):
            bases_a = range(int(b[5]),int(b[4]))
        if int(b[6]) < int(b[7]):
            bases_b = range(int(b[6]),int(b[7]))
        if int(b[6]) > int(b[7]):
            bases_b = range(int(b[7]),int(b[6]))
        a_s = set(bases_a)
        output.write(b[0]+'\t'+b[2]+'\t'+str(float(len(a_s.intersection(bases_b))/float(len(bases_a))))+'\n')
                
