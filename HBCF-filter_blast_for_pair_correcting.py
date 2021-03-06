#!/bin/usr/env python 
import sys
import argparse

parser = argparse.ArgumentParser(description='''Export allele and ohnologues tables indicating the unpaired gene and the new allele/ohnologue with scaffold positional information etc''')
parser.add_argument('BLAST_infile', help = 'the name of the BLAST file - must be outfmt 10 as outlined in the HBCF wiki')
parser.add_argument('-unpaired', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/UNPAIRED_GENES.txt',  help = 'a list of names corresponding to the id of the unpaired gene created by HBCF-gene_and_pair_mining.py')
parser.add_argument('-genes', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', help = 'The full path to the genes.all file of the A.vaga genome v2 or a tab separated table: first column [scaffold id]; second column [gene_id]; third column [gene start position]; fourth column [gene end position] e.g. [av1][GSADVT00000034001][123][225]')
parser.add_argument('-pairs', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.all', help = 'The full path to the pairs.all file of the A.vaga genome v2 or a tab separated table : first column [random]; second column[random]; third column [gene_1]; fourth column [gene_2] e.g. [0][0][GSADVT00000227001][GSADVT00022299001]')
parser.add_argument('-ohno_ident', default = '70.0', type = float, help = 'The minimum percent identitiy to be considered for ohnologue assignment')
parser.add_argument('-allelic_ident', default = '90.0', type = float, help = 'The minimum percent identity to be considered for alleleic assignment')
parser.add_argument('-length', default = '0.75', type = float, help = 'The minimum percentage length that the hit coveres the query sequence')
parser.add_argument('-o', help = 'the directory to write output files')
args = parser.parse_args()

infile = open(args.BLAST_infile, 'rU')
genes =open(args.genes, 'rU')
pairs =open(args.pairs, 'rU')
unpaired =open(args.unpaired,'rU')


blast =[]
for line in infile: # parse the blast file
    x= line.strip().split(",")
    header = x[0].split("_")
    blast.append([header[0], header[1]+'_'+header[2], header[3].strip("bp"), header[4], header[5], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]])
    
hits_to_self = []
allele_hits = []
ohno_hits = []

for line in blast:

    if int(line[7]) < 250:
        next
    elif line[1] == line[5] and int(line[3]) == int(line[12])-1 and float(line[6]) == 100.0:
        hits_to_self.append(line)
    elif float(float(line[7])/float(line[2])) > args.length and line[1] != line[5] and float(float(line[6])) < args.allelic_ident and float(float(line[6])) > args.ohno_ident:
        ohno_hits.append(line)
    elif float(float(line[7])/float(line[2])) > args.length and line[1] == line[5] and int(line[3]) != int(line[12])-1 and float(float(line[6])) < args.allelic_ident and float(float(line[6])) > args.ohno_ident:
        ohno_hits.append(line)
    elif float(float(line[7])/float(line[2])) > args.length and line[1] == line[5] and int(line[3]) != int(line[12])-1 and float(float(line[6])) > args.allelic_ident:
        allele_hits.append(line)
    elif float(float(line[7])/float(line[2])) > args.length and line[1] != line[5] and float(float(line[6])) > args.allelic_ident:
        allele_hits.append(line)

print("whew, you recovered %d alleles" % (len(allele_hits)))
print("awesome, you recovered %d ohnos" % (len(ohno_hits)))
print("hmmm, you found %d hits to self" % (len(hits_to_self)))

gene_dict = {} # Create a dictionary for the genes in the A.vaga consortium data.
for gene in genes:
    g = gene.strip().split("\t")
    scaff = g[0].strip("av")
    new_scaff = 'scaffold_'+scaff
    gene_dict[g[1]] = [new_scaff, g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]]

    
def check_for_gene_overlap(unpaired_list, dict): # This function will search through gene_dict for a matching gene regrion to the quality_hits line
    d = 0
    for key in dict.keys():
        if dict[key][0] == unpaired_list[5]: # Look for the matching scaffold to the line
            a = range(int(unpaired_list[12]), int(unpaired_list[13])) 
            b = range(int(dict[key][3]), int(dict[key][4]))
            a_s = set(a)
            new_gene = [] # this is a list holder that we can store relevant information if we find a match to an existing gene.
            if len(a_s.intersection(b)) > 1: # Im starting with one here, but I don't kow what constitutes an already existing gene. Probably more than the ovelap of one position.
                c = dict[key][2].split("T")
                if len(str(int(c[1])+1)) == 8:
                    n = "000"
                if len(str(int(c[1])+1)) == 7:
                    n = "0000"
                if len(str(int(c[1])+1)) == 6:
                    n = "00000"
                if len(str(int(c[1])+1)) == 5:
                    n = "000000"
                if len(str(int(c[1])+1)) == 4:
                    n = "0000000"
                new_name = str(d)+"_"+c[0]+'T'+n+str(int(c[1])+1)
                print new_name
                               # the new gene contains the following variables [0:Unpaired_gene], [1:scaffold], [2:new_gene_name], [3:existing_gene_scaffold], [4:existing_gene_start], [5:existing_gene_end], [4:query_gene_start], [5:query_gene_stop], [6:percent_id], [7:length], [8:evalue]
                new_gene.append([unpaired_list[0], unpaired_list[1], new_name, dict[key][0], dict[key][3], dict[key][4], unpaired_list[12], unpaired_list[13], unpaired_list[6], unpaired_list[7], unpaired_list[14]])
                return new_gene
        d += 1

####################################

allele_new_genes = []                 
allele_existing_genes = []
d = 0
for p in allele_hits:
    w = check_for_gene_overlap(p, gene_dict)
    if w != None:
        allele_existing_genes.append(w[0])
    if w == None: # These are the genes that have no overlap with existing genes.
        c = p[0].split("T")
        if len(str(int(c[1])+1)) == 8:
             n = "000"
        if len(str(int(c[1])+1)) == 7:
            n = "0000"
        if len(str(int(c[1])+1)) == 6:
            n = "00000"
        if len(str(int(c[1])+1)) == 5:
            n = "000000"
        if len(str(int(c[1])+1)) == 4:
            n = "0000000"
        new_name = str(d)+"_"+c[0]+'T'+n+str(int(c[1])+500)
        print new_name
        # appending a list of [0:Unpaired_gene], [1:Unpaired_scaffold], [2:new_gene_name], [3: hit_scaffold], [4: hit_start], [5: hit_stop], [6: percent_ID], [7:length], [8: evalue] 
        allele_new_genes.append([p[0], p[1], new_name, p[5], p[12], p[13], p[6], p[7], p[14]])
    d += 1

#unpaired_genes = open('UNPAIRED_GENE_COORDINATES.txt', 'w')
allele_ng_outfile = open(args.o+'/ALLELES_NEW_GENES.txt', 'w')
allele_eg_outfile = open(args.o+'/ALLELES_EXISTING_GENES.txt', 'w')
allele_ng_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"New_Gene_Scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evlaue"+"\n")
allele_eg_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"Existing_gene_scaffold"+"\t"+"NEW_Gene_Start"+"\t"+"NEW_Gene_End"+"\t"+"Existing_gene_start"+"\t"+"Existing_gene_end"+"\t"+"percent_id"+"\n")
#unpaired_genes.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"Gene_Start"+"\t"+"Gene_End"+"\n")

#for line in hits_to_self:
#    unpaired_genes.write(line[0]+"\t"+line[1]+"\t"+line[10]+"\t"+line[11]+"\n")
new_gene_names = []
for line in allele_new_genes:
    allele_ng_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
n_count = 0
for line in allele_existing_genes:
    allele_eg_outfile.write(line[0]+"\t"+line[1]+"\t"+str(n_count)+'-'+line[2]+"\t"+line[3]+"\t"+line[6]+"\t"+line[7]+"\t"+line[4]+"\t"+line[5]+"\t"+line[8]+"\n")
    n_count += 1

############################ 

ohno_new_genes = []
ohno_existing_genes = []
d = 0
for p in ohno_hits:
    w = check_for_gene_overlap(p, gene_dict)
    if w != None:
        ohno_existing_genes.append(w[0])
    if w == None: # These are the genes that have no overlap with existing genes.
        c = p[0].split("T")
        if len(str(int(c[1])+1)) == 8:
             n = "000"
        if len(str(int(c[1])+1)) == 7:
            n = "0000"
        if len(str(int(c[1])+1)) == 6:
            n = "00000"
        if len(str(int(c[1])+1)) == 5:
            n = "000000"
        if len(str(int(c[1])+1)) == 4:
            n = "0000000"
        new_name = str(d)+"_"+c[0]+'T'+n+str(int(c[1])+500)
        print new_name
        ohno_new_genes.append([p[0],p[1], new_name, p[5], p[12], p[13], p[6], p[7],p[14]])
    d += 1


#unpaired_genes = open('UNPAIRED_GENE_COORDINATES.txt', 'w')
ohno_ng_outfile = open(args.o+'/OHNO_NEW_GENES.txt', 'w')
ohno_eg_outfile = open(args.o+'/OHNO_EXISTING_GENES.txt', 'w')
ohno_ng_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"New_Gene_Scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evlaue"+"\n")
ohno_eg_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"Existing_gene_scaffold"+"\t"+"Exisiting_Gene_Start"+"\t"+"Existing_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evalue"+"\n")
#unpaired_genes.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"Gene_Start"+"\t"+"Gene_End"+"\n") 


#for line in hits_to_self:
#    unpaired_genes.write(line[0]+"\t"+line[1]+"\t"+line[10]+"\t"+line[11]+"\n")
for line in ohno_new_genes:
    ohno_ng_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
o_counts = 0
for line in ohno_existing_genes:
    ohno_eg_outfile.write(line[0]+"\t"+line[1]+"\t"+str(o_counts)+'-'+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
    o_counts += 1
