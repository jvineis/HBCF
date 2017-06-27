#!/bin/usr/env python 
import sys
import argparse

parser = argparse.ArgumentParser(description='''Export allele and ohnologues tables indicating the unpaired gene and the new allele/ohnologue with scaffold positional information etc''')
parser.add_argument('BLAST_infile', help = 'the name of the BLAST file - must be outfmt 10 as outlined in the HBCF wiki')
parser.add_argument('--genes', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all', help = 'The full path to the genes.all file of the A.vaga genome v2 or a tab separated table: first column [scaffold id]; second column [gene_id]; third column [gene start position]; fourth column [gene end position] e.g. [av1][GSADVT00000034001][123][225]')
parser.add_argument('--pairs', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.all', help = 'The full path to the pairs.all file of the A.vaga genome v2 or a tab separated table : first column [random]; second column[random]; third column [gene_1]; fourth column [gene_2] e.g. [0][0][GSADVT00000227001][GSADVT00022299001]')
args = parser.parse_args()

infile = open(args.BLAST_infile, 'rU')
genes =open(args.genes, 'rU')
pairs =open(args.pairs, 'rU')

blast =[]
for line in infile: # parse the blast file
    x= line.strip().split(",")
    header = x[0].split("_")
    blast.append([header[0],header[1]+'_'+header[2], header[3].strip("bp"),  x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8],x[9],x[10],x[11]])
    
hits_to_self = []
allele_hits = []
ohno_hits = []
for line in blast:
    if float(line[5]) < 300.0:
        next
    elif float(float(line[5])/float(line[2])) > 0.9 and line[1] != line[3] and line[2] != line[5] and float(float(line[4])) > 93.0:
        allele_hits.append(line)
    elif float(float(line[5])/float(line[2])) > 0.9 and line[1] != line[3] and line[2] != line[5] and float(float(line[4])) < 93.0 and float(float(line[4])) > 75.0:
        ohno_hits.append(line)
    elif float(float(line[5])/float(line[2])) == 1.0 and line[1] == line[3] and line[2] == line[5]:
        hits_to_self.append(line)

#for line in allele_hits:
#    print line

gene_dict = {}
for gene in genes:
    g = gene.strip().split("\t")
    scaff = g[0].strip("av")
    new_scaff = 'scaffold_'+scaff
    gene_dict[g[1]] = [new_scaff, g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]]
    
def check_for_gene_overlap(unpaired_list, dict): # This function will search through gene_dict for a matching gene regrion to the quality_hits line
    for key in dict.keys():
        if dict[key][0] == unpaired_list[3]: # Look for the matching scaffold to the line
#            print unpaired_list
            a = range(int(unpaired_list[10]), int(unpaired_list[11])) 
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
                new_name = c[0]+'T'+n+str(int(c[1])+1)
                new_gene.append([unpaired_list[0], unpaired_list[1],new_name,dict[key][0],unpaired_list[10],unpaired_list[11],dict[key][3],dict[key][4], unpaired_list[4], unpaired_list[5], unpaired_list[12]])
                print "we found gene %s scaffold %s with a start %d and stop %d that matches %s scaffold with a start %s and end %s" %(gene_dict[key][2],gene_dict[key][0],int(gene_dict[key][3]),int(gene_dict[key][4]),unpaired_list[3], unpaired_list[10], unpaired_list[11]), unpaired_list[4],unpaired_list[5],unpaired_list[12]
                return new_gene


####################################

allele_new_genes = []                 
allele_existing_genes = []
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
        new_name = c[0]+'T'+n+str(int(c[1])+500)
        allele_new_genes.append([p[0],p[1], new_name, p[3], p[10], p[11], p[4], p[5], p[12]])


#unpaired_genes = open('UNPAIRED_GENE_COORDINATES.txt', 'w')
allele_ng_outfile = open('ALLELES_NEW_GENES.txt', 'w')
allele_eg_outfile = open('ALLELES_EXISTING_GENES.txt', 'w')
allele_ng_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"New_Gene_Scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evlaue"+"\n")
allele_eg_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"Existing_gene_scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evalue"+"\n")
#unpaired_genes.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"Gene_Start"+"\t"+"Gene_End"+"\n")

#for line in hits_to_self:
#    unpaired_genes.write(line[0]+"\t"+line[1]+"\t"+line[10]+"\t"+line[11]+"\n")
for line in allele_new_genes:
    allele_ng_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
for line in allele_existing_genes:
    allele_eg_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")


############################ 

ohno_new_genes = []
ohno_existing_genes = []
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
        new_name = c[0]+'T'+n+str(int(c[1])+500)
        ohno_new_genes.append([p[0],p[1], new_name, p[3], p[10], p[11], p[4], p[5],p[12]])


#unpaired_genes = open('UNPAIRED_GENE_COORDINATES.txt', 'w')
ohno_ng_outfile = open('OHNO_NEW_GENES.txt', 'w')
ohno_eg_outfile = open('OHNO_EXISTING_GENES.txt', 'w')
ohno_ng_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"New_Gene_Scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"\t"+"evlaue"+"\n")
ohno_eg_outfile.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"New_Gene_Name"+"\t"+"Existing_gene_scaffold"+"\t"+"New_Gene_Start"+"\t"+"New_Gene_End"+"\t"+"percent_ID"+"\t"+"length"+"t"+"evalue"+"\n")
#unpaired_genes.write("Unpaired_gene"+"\t"+"Unpaired_scaffold"+"\t"+"Gene_Start"+"\t"+"Gene_End"+"\n") 



#for line in hits_to_self:
#    unpaired_genes.write(line[0]+"\t"+line[1]+"\t"+line[10]+"\t"+line[11]+"\n")
for line in ohno_new_genes:
    ohno_ng_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
for line in ohno_existing_genes:
    ohno_eg_outfile.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n")
