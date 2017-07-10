#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Extract the cds of a genewise output and the allelic cds from the A.vaga genome and write both sequences to a single fasta file for alignment with muscle''')
parser.add_argument('fasta_file_name', help = 'the name of the fasta file containing the sequence hit by the unpaired gene')
parser.add_argument('--genewise', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/GENEWISE/', help='The directory containing the genewise output.  Files within the directory must have names that link it to the fasta file input name e.g. fasta_file.genewise.out.txt')
parser.add_argument('--o_pep', default='.', help='The location for output of the peptide fasta')
parser.add_argument('--o_nuc', default='.', help='The location for output of the nucleotide fasta')
parser.add_argument('--ref_fasta', default='/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.cds.fa', help = 'The location of the fasta file containing reference sequences to align with the cds produced by genewise. Should contain the cds of the original unpaired gene')
parser.add_argument('--ref_pep', default='/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.pep.fa', help = 'The location of the peptide fasta file containing reference sequences to align with the peptide produced by genewise.  Should contain the peptide of the original unpaired gene')
args = parser.parse_args()

genewise = args.fasta_file_name # name of the scaffold fasta file containing sequence hit by unpaired genes
genewise_in = open(args.genewise+genewise+'.genewise.out.txt', 'rU') # open the genewise output
outfile_fasta = open(args.o_nuc+genewise+'.cds_to_align.fa', 'w') # open the outfile to write the cds sequences
outfile_pep = open(args.o_pep+genewise+'.pep_to_align.fa', 'w') # open the outfile to wirte the peptide sequences
file_name = genewise.split("_") # separating the characters of the fasta file neede to select from the fasta containg the cds
#print file_name
gsadv = file_name[0][0:14]+'001' # the characters needed to select the cds from the Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa
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

k = [i for i, j in enumerate(genewise_contents) if '>' in j] # record the positions in the file where the peptide and fasta sequences are
peps_and_fasta = {} # create a dictionary to record the sequence and name

print("there are %s peptide and nucleotide seqs, and the file is %s lines long" %(len(k), len(genewise_contents))) # just so we know that there are lines in the file with data

line_count = len(genewise_contents) # record the number representing the last line of the file
n = 0 # a counter
for start in k: # for each of the positions contaning a sequence name build a dictionary for the name of the sequence and concatenate the lines of sequence data into a single sequence for each entry
    if n < int(len(k)-1): # this part of the loop allows us to stay in bounds of the number of start positions in "k" 
        name =  genewise_contents[k[n]]
        seq = "".join(genewise_contents[int(k[n])+1:int(k[n+1])])
        peps_and_fasta[n] = [name,seq]
        n += 1
    else:
        name = genewise_contents[k[len(k)-1]]
        seq = "".join(genewise_contents[int(k[len(k)-1])+1:int(line_count)])
        peps_and_fasta[len(k)-1] = [name,seq]


dict_keys = [] # a list of the keys in the dictionary created above
for key in peps_and_fasta.keys():
    dict_keys.append(key)

fasta_files = {}
peptide_files = {}
bases = ["A","G","C","T"]

for key in peps_and_fasta.keys(): # This bit of code below takes care of the instances where the nucleotide sequence is divided into more than one part by an intron. I'm taking the two parts and pasting them together into a single file 
    if "sp" in peps_and_fasta[key][0]:
        seq = peps_and_fasta[key][1].split("//")
        fasta_files[key] = [peps_and_fasta[key][0],len(seq[0]), seq[0]]

for key in peps_and_fasta.keys(): # Here we write the peptide sequences to a dictionary
     if "pep" in peps_and_fasta[key][0]:
        pep = peps_and_fasta[key][1].split("//")
        peptide_files[key] = [peps_and_fasta[key][0], len(pep[0]), pep[0]]

peptide_seqs = [] # Put the peptide sequences from the peptide dictionary into a single list
for key in peptide_files:
    peptide_seqs.append(peptide_files[key][2])

nucleotide_seqs = [] # Put the nucleotide sequences from the peptide dictionary into a single list
for key in fasta_files.keys():
    nucleotide_seqs.append(fasta_files[key][2])

big_fasta_seq = "".join(nucleotide_seqs)
print("your seq list", nucleotide_seqs)
big_pep_seq = "".join(peptide_seqs)
print("your nucleotide seq", big_fasta_seq)
print("your peptide seq", big_pep_seq)
genewise_fasta_header = file_name[0][0:14]+'_genewise'

outfile_fasta.write(">"+genewise_fasta_header+"\n"+big_fasta_seq)
outfile_pep.write(">"+genewise_fasta_header+"_pep"+"\n"+big_pep_seq)



