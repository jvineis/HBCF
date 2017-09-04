#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Extract the cds of a genewise output and the allelic cds from the A.vaga genome and write both sequences to a single fasta file for alignment with muscle''')
parser.add_argument('fasta_file_name', help = 'the name of the fasta file containing the sequence hit by the unpaired gene')
parser.add_argument('matching_locations', help = 'the table containg the unpaired hits to new or existing genes - ususally called ALLELES_NEW_GENES.txt or ALLELES_EXISTING_GENES.txt')
parser.add_argument('--genewise', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/GENEWISE/', help='The directory containing the genewise output.  Files within the directory must have names that link it to the fasta file input name e.g. fasta_file.genewise.out.txt')
parser.add_argument('--o_pep', default='.', help='The location for output of the peptide fasta')
parser.add_argument('--o_nuc', default='.', help='The location for output of the nucleotide fasta')
parser.add_argument('--ref_fasta', default='/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.cds.fa', help = 'The location of the fasta file containing reference sequences to align with the cds produced by genewise. Should contain the cds of the original unpaired gene')
parser.add_argument('--ref_pep', default='/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.pep.fa', help = 'The location of the peptide fasta file containing reference sequences to align with the peptide produced by genewise.  Should contain the peptide of the original unpaired gene')
parser.add_argument('--score', default=20.0, type=float, help='the minimum score to accept the peptide sequence')
args = parser.parse_args()

matching_locations = open(args.matching_locations, 'rU')
genewise = args.fasta_file_name # name of the scaffold fasta file containing sequence hit by unpaired genes
genewise_in = open(args.genewise+genewise+'.genewise.out.txt', 'rU') # open the genewise output
outfile_fasta = open(args.o_nuc+genewise+'.cds_to_align.fa', 'w') # open the outfile to write the cds sequences
outfile_pep = open(args.o_pep+genewise+'.pep_to_align.fa', 'w') # open the outfile to wirte the peptide sequences
file_name = genewise.split("_") # separating the characters of the fasta file neede to select from the fasta containg the cds
gsadv_f = file_name[1][0:14]+'001' # the characters needed to select the cds from the Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa
for line in matching_locations:
    x = line.strip().split("\t")
    target = line.strip().split("\t")[2].split("_")[1][0:14]+'001'
    if gsadv_f == target:
        gsadv = x[0]
        print gsadv, target
        
print("gene analyzed",gsadv)
cds_fasta = open(args.ref_fasta, 'rU')
pep_fasta = open(args.ref_pep, 'rU')


with cds_fasta as f: # Parse the fasta file corresponding to the cds found in the consortium files
    cds = SeqIO.parse(f, "fasta")
    for rec in cds:
        if gsadv in rec.id:
            outfile_fasta.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n') #Write the header and sequence to the outfile

with pep_fasta as f:
    pep = SeqIO.parse(f, "fasta")
    for rec in pep:
        if gsadv in rec.id:
            outfile_pep.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n') # Write the pep header and sequence to the outfile


genewise_contents = []# The rest of the code below is designed to pull the longest sequence out of the genewise output and add it to the existing outfile
for line in genewise_in: # Parse the genewise infile into individual lines
    x = line.strip("\n")
    genewise_contents.append(x)

a = [i for i, j in enumerate(genewise_contents) if 'Score' and 'alignment' in j] # record the position in the file where the Score is reported 
genewise_contents_score_filtered = [] # Here we are making a subset of the "pretty" genewise output based on whether or not the score of the peptide hit is above a certain threshold (usually 20).  This is empty list to catch the results
count = 0  # a counter
for i in a: # for the positions of the genewise output that contains the score information e.g. [2,36,55]
    info = genewise_contents[i].split(" ") # split the line by spaces
    if count < int(len(a)-1) and float(info[1]) > float(args.score): # if 
        genewise_contents_score_filtered.append(genewise_contents[i:a[int(count)+1]])
    elif float(info[1]) > float(args.score):
        genewise_contents_score_filtered.append(genewise_contents[i:len(genewise_contents)])
    count += 1

#print genewise_contents_score_filtered

def dna_and_pep_puller(genewise_fraction_to_analyze):
    dna_and_pep = {}
    k = [i for i, j in enumerate(genewise_fraction_to_analyze) if 'GSADVT' in j]
    line_count = len(genewise_fraction_to_analyze)
    n = 0
    for start in k:
        if n < int(len(k)-1):
            name = genewise_fraction_to_analyze[k[n]]
            seq = "".join(genewise_fraction_to_analyze[int(k[n])+1:int(k[n+1])])
            dna_and_pep[n] = [name,seq]
            n += 1
        else:
            name = genewise_fraction_to_analyze[k[len(k)-1]]
            seq = "".join(genewise_fraction_to_analyze[int(k[len(k)-1])+1:int(line_count)])
            dna_and_pep[len(k)-1] = [name,seq]
    return dna_and_pep


dna_and_pep_master = []
for i in genewise_contents_score_filtered: # This loops through each of the individually reported DNA and PEP hits from genewise and pulls out the sequence information (DNA and PEP) from the rest of the stuff
    d = dna_and_pep_puller(i) # Run the dna_and_pep_puller function on the report which returns the name of the seq and the sequence
    for key in d.keys(): # for each of the sequence lines   
        dna_and_pep_master.append(d[key][0:len(d[key])])


fasta_files = {}
peptide_files = {}
bases = ["A","G","C","T"]

line = 0
for dna in dna_and_pep_master: # This bit of code below takes care of the instances where the nucleotide sequence is divided into more than one part by an intron. I'm taking the two parts and pasting them together into a single file 
    if "sp" in dna[0]:
        seq = dna[1].split("//")
        fasta_files[line] = [dna[0],len(seq[0]), seq[0]]
    line += 1

line = 0
for pep in dna_and_pep_master: # Here we write the peptide sequences to a dictionary
    if "pep" in pep[0]:
        pep = pep[1].split("//")
        peptide_files[line] = [dna_and_pep_master[0], len(pep[0]), pep[0]]
    line += 1

peptide_seqs = [] # Put the peptide sequences from the peptide dictionary into a single list
for key in peptide_files:
    peptide_seqs.append(peptide_files[key][2])

nucleotide_seqs = [] # Put the nucleotide sequences from the nucleotide dictionary into a single list
for key in fasta_files.keys():
    nucleotide_seqs.append(fasta_files[key][2])

big_fasta_seq = "".join(nucleotide_seqs)
big_pep_seq = "".join(peptide_seqs)
print("your nucleotide seq", big_fasta_seq)
print("your peptide seq", big_pep_seq)
genewise_fasta_header = file_name[0]+"_"+file_name[1]

outfile_fasta.write(">"+genewise_fasta_header+"\n"+big_fasta_seq)
outfile_pep.write(">"+genewise_fasta_header+"_pep"+"\n"+big_pep_seq)



