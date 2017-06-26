#!/usr/bin/env python
import sys
import os
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import csv
import StringIO
from glob import glob
###Sequence fasta file - this should be the same fasta file that you used to create the separate allele fasta files for 
### mapping (allele_pairs_one.fa and allele_pairs_two.fa) using "split_pairs_to_fa.py".
handle = open(sys.argv[1], 'rU')
refseqs = SeqIO.parse(handle, "fasta")
###Pairs file - this file contains the ID pairs in tab delimited form and is the same file used in the "split_pairs_to_fa.py"
## command.
pairs = open(sys.argv[2], 'rU')


## Make a single list with two columns parsed
seqlist = []
for line in pairs:
    x = line.strip("\n").split("\t")
    seqlist.append([x[0],x[1]])

## We need to make a dictionary for the fasta file -  sequence id as the key and each SeqIO record as the value
seq_dict = {}
for rec in refseqs:
    seq_dict[rec.id] = rec

## Create a list for the first column of the paired sequences list
one_list = []
for line in seqlist:
    one_list.append(line[0])
print "%d seqids in list one" % len(one_list)
## Create another list for the second column of the paired sequences list
two_list = []
for line in seqlist:
    two_list.append(line[1])
print "%d seqids in list two" % len(two_list)
## Search through the dictionary for the ids in the first list and append the hits to the "one_recs" list
one_recs = []
for line in one_list:
    if seq_dict.has_key(line):
        one_recs.append(seq_dict[line])
print "%d recs in seqrec one" % len(one_recs)
#Search through the dictionary for the ids in the second list and append the hits to the "two_recs" list
two_recs = []
for line in two_list:
    if seq_dict.has_key(line):
        two_recs.append(seq_dict[line])
print "%d recs in seqrec two" % len(two_recs)
## Merge the two lists of seq records together
z = 1
merged_recs = zip(one_recs,two_recs)
print "%d recs in merged_recs" % len(merged_recs)
for line in merged_recs:
    outfile = line[0]
    outfile = outfile.id+"_"+str(z)
    outfile = outfile+"-pairs.fa"
    seq_handle = open(outfile,'w')
    SeqIO.write([line[0],line[1]],seq_handle, "fasta")
    z+=1
seq_handle.close()

## This section of code creates the muscle command line and runs the command
## for each alignment.

x=1
for filename in glob('*pairs.fa'):
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
matching_pairs_snps_list = []

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
        ## 3. create a variable for the pair1 sequence list
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
                    pair_one_snps_to_search.append([a[0].id,seq_dict_one[pos],pair1_sequence[pos], pair2_sequence[pos]])
                    pair_two_snps_to_search.append([a[1].id,seq_dict_two[pos],pair1_sequence[pos], pair2_sequence[pos]])
                    matching_pairs_snp_list.append([[str(a[0].id)+"_"+str(seq_dict_one[pos])+"_"+str(pair1_sequence[pos])+"_"+str(pair2_sequence[pos])],[str(a[1].id)+"_"+str(seq_dict_two[pos])+"_"+str(pair1_sequence[pos])+"_"+str(pair2_sequence[pos])]])
    print "%d Pairs analyzed" % (y)
    y += 1

with outfile_one as f:
    w = csv.writer(f, dialect = 'excel-tab')
    w.writerows(pair_one_snps_to_search)

with outfile_two as f:
    w = csv.writer(f, dialect = 'excel-tab')
    w.writerows(pair_two_snps_to_search)

with outfile_three as f:
    w = csv.writer(f, dialect = 'excel-tab')
    w.writerows(matching_pairs_snp_list)

outfile_one.close()
outfile_two.close()
pairs.close()
handle.close()
