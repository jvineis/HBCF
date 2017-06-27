#!/usr/bin/env python
import sys
import os
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import csv
import StringIO
from glob import glob
print("This script is making awesome outputs for the allele_name_1 snp_position_1 base_1 allele_name_2 snp_position_2 base_2.  Make sure that you are running it from the location where the fasta files exist or nothing great will happen")

x=1
for filename in glob('*.fa'):
    with open(filename) as f:
        output = str(filename)
        output += '-aligned.aln'
        in_file = str(filename)
        muscle_cline = MuscleCommandline(input=in_file, out=output, clw=True)
#        print(muscle_cline)
        stdout, stderr = muscle_cline()
        print "%d Pairs aligned" % (x)
        x += 1

## Here we scan through each alignment and find the locations in the alignment
## where there is a mismatch.  Gaps in the alignment are ignored.

## An outfile for the contig name, snp position, and bases for pair one
outfile_one = open('PAIRS_ONE_SNPS.txt', 'w')
## An outfile for the contig name, snp position, and bases for pair two
outfile_two = open('PAIRS_TWO_SNPS.txt', 'w')
## An outfile for the matching pairs of snps
outfile_three = open('MATCHING_PAIRS_SNPS.txt', 'w')

## A list to hold the values of the snp position search for pair one
pair_one_snps_to_search = []
## A list to hold the values of the snp position search for pair two
pair_two_snps_to_search = []
## A list to hold the matiching values of the snp pairs
matching_pairs_snp_list = []

y = 1
n_one = 1
n_two = 1
for filename in glob('*.aln'):
    with open(filename) as a:
        in_file = str(filename)
        a = AlignIO.read(in_file, "clustal")
        sequences = ["A","C","G","T"]
        gaps = ["-"]
        ## Here we need a variable to keep track of the position in the alignment
        ## And we one that we can match with the reference position needed to
        ## search in the variant tables produced by CLC varaint caller

        ## lets try creating a dictionary for the alignment pairs that
        ## uses the position in the alignment as the key and the position
        ## in the reference as the value. The difference in position is due to the
        ## gap "-" character.  These do not exist in the reference

        ## 1. create a dictionary variable for alignment 1
        seq_dict_one = {}
        ## 2. create a dictionary variable for alignment 2
        seq_dict_two = {}
        pair1_sequence = a[0].seq
        ## 4. create a variable for the pari2 sequence list
        pair2_sequence = a[1].seq
        n_one = 0
        for pos in range(len(a[0].seq)):
            if pair1_sequence[pos] in sequences:
                n_one += 1
                seq_dict_one[pos] = n_one
            else:
                seq_dict_one[pos] = pos

#        print seq_dict_one
        n_two = 0
        for pos in range(len(a[0].seq)):
            if pair2_sequence[pos] in sequences:
                n_two +=1
                seq_dict_two[pos] = n_two
            else:
                seq_dict_two[pos] = pos
#        print seq_dict_two

        for pos in range(len(a[0].seq)):
            if pair1_sequence[pos] in sequences and pair2_sequence[pos] in sequences:
                if pair1_sequence[pos] != pair2_sequence[pos]:
                    name_1 = a[0].id.split("_")[0]
                    name_2 = a[1].id.split(":")[0]
                    pair_one_snps_to_search.append([name_1, str(seq_dict_one[pos]), str(pair1_sequence[pos]), str(pair2_sequence[pos])])
                    pair_two_snps_to_search.append([name_1, str(seq_dict_two[pos]), str(pair1_sequence[pos]), str(pair2_sequence[pos])])
                    matching_pairs_snp_list.append([name_1, str(seq_dict_one[pos]), str(pair1_sequence[pos]), str(pair2_sequence[pos]), name_2, str(seq_dict_two[pos]), str(pair1_sequence[pos]), str(pair2_sequence[pos])])

    print "%d Pairs analyzed" % (y)
    y += 1

for line in pair_one_snps_to_search:
    outfile_one.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n")
for line in pair_two_snps_to_search:
    outfile_two.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n")
for line in matching_pairs_snp_list:
    outfile_three.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\n")
