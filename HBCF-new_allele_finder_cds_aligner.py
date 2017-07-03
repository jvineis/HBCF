#!/usr/bin/env python

import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Extract the cds of a genewise output and the allelic cds from the A.vaga genome and write both sequences to a single fasta file for alignment with muscle''')
parser.add_argument('fasta_file', help = 'the fasta file containing the sequence hit by the unpaired gene')
parser.add_argument('--genewise', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/GENEWISE/', help='The directory containing the genewise output.  Files within the directory must have names that link it to the fasta file input name e.g. fasta_file.genewise.out.txt')
parser.add_argument('--peptides', default='/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/NEW_ALLELES/CDS_ALIGNMENT/', help='The location for output of the fasta files')
parser.add_argument('--ref_fasta', default='/groups/rotifer/Avgenome/Genoscope/v2/Adineta_vaga_v2.0.annot.cds.fa', help = 'The location of the fasta file containing reference sequences to align with the cds produced by genewise. Should contain the cds of the original unpaired gene')
args = parser.parse_args()

genewise = args.fasta_file # name of the scaffold fasta file containing sequence hit by unpaired genes
genewise_in = open(args.genewise+genewise+'.genewise.out.txt', 'rU') # open the genewise output
outfile = open(args.peptides+genewise+'.cds_to_align.fa', 'w') # open the outfile to write the cds sequences
file_name = genewise.split("_") # separating the characters of the fasta file neede to select from the fasta containg the cds
#print file_name
gsadv = file_name[0][0:14]+'001' # the characters needed to select the cds from the Adineta_vaga_v2.0.annot.transcripts.forAnvio.fa
print(gsadv)
cds_fasta = open(args.ref_fasta, 'rU')

with cds_fasta as f: # Parse the fasta file corresponding to the cds found in the consortium files
    cds = SeqIO.parse(f, "fasta")
    for rec in cds:
        if gsadv in rec.id:
            outfile.write('>'+str(rec.id)+'\n'+str(rec.seq)+'\n') #Write the header and sequence to the outfile

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

for key in peps_and_fasta.keys(): # This bit of code below takes care of the instances where the nucleotide sequence is divided into more than one part by an intron. I'm taking the two parts and pasting them together into a single file 
    if "sp" in peps_and_fasta[key][0] and "//" in peps_and_fasta[key][1]:
        seq = peps_and_fasta[key][1].strip("//")
        fasta_files[key] = [peps_and_fasta[key][0],len(seq), seq]

    elif "sp" in peps_and_fasta[key][0] and "//" not in peps_and_fasta[key][1]:
        fasta_seq = []
        fasta_seq.append(peps_and_fasta[key][1])
        for i in dict_keys:
            if int(i) > int(key) and "//" not in peps_and_fasta[i][1]:
                fasta_seq.append(peps_and_fasta[i][1])
            elif "//" in peps_and_fasta[i][1] and "sp" in peps_and_fasta[i][0]:
                fasta_seq.append(peps_and_fasta[i][1].strip("//"))
                break
        fasta_files[key] = [peps_and_fasta[key][0], len("".join(fasta_seq)), "".join(fasta_seq)]

for key in peps_and_fasta.keys(): # Here we write the peptide sequences to a dictionary althought this code doesn't do anything with this at the moment.
     if "pep" in peps_and_fasta[key][0]:
        pep = peps_and_fasta[key][1].strip("//")
        peptide_files[key] = [peps_and_fasta[key][0], len(pep), pep]
        
#print(fasta_files)
#print(peptide_files)

nucleotide_seqs = []
for key in fasta_files.keys():
    nucleotide_seqs.append(fasta_files[key][2])

big_seq = max(nucleotide_seqs, key = len)
genewise_fasta_header = file_name[0][0:14]+'_genewise'

outfile.write(">"+genewise_fasta_header+"\n"+big_seq)



