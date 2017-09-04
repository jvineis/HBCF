#!/usr/bin/env python
# The purpose of this script is to create a bash script that can be run to examine the presence of snps in all splits of a gene for an individual scaffold.
import sys

scaf_num = sys.argv[1]

outfile_name = 'run_assign_var_pos_'+scaf_num+'.shx'
outfile = open(outfile_name, 'w')
outfile.write("#!/bin/bash\n")
outfile.write("for i in `cat av%s_GENE_list.txt`; do python ~/scripts/HBCF-anvio_assign_variant_positions.py scaffold_%s_genes_per_split.txt $i; done\npython ~/scripts/HBCF-anvio_assign_variant_positions.py scaffold_%s intron\ncat *snps_per_split.txt intron_snps_per_split.txt | sort > Scaffold_%s_introns_and_snps.txt" %(scaf_num, scaf_num, scaf_num,scaf_num))
              
#for i in `cat av1_GENE_list.txt`; do python ~/scripts/HBCF-anvio_assign_variant_positions.py scaffold_1_genes_per_split.txt $i; done
#python ~/scripts/HBCF-anvio_assign_variant_positions.py scaffold_1_genes_per_split.txt intron
#cat *snps_per_split.txt intron_snps_per_split.txt | sort > Scaffold_1_introns_and_snps.txt
