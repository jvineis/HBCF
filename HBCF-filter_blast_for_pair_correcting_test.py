#!/bin/usr/env python
import sys

infile = open(sys.argv[1], 'rU')
genes =open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/genes.all')
pairs =open('/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.all')

blast =[]
for line in infile: # parse the blast file
    x= line.strip().split(",")
    header = x[0].split("_")
    blast.append([header[0],header[1]+'_'+header[2], header[3].strip("bp"),  x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8],x[9],x[10],x[11]])

hits_to_self = []
#quality_hits = [] # identify quality hits that are at least 90 percent of lenght and are not a hit to its self
allele_hits = []
ohno_hits = []
for line in blast:
    if float(line[5]) < 300.0:
#        print("too short", line)
        next 
    elif float(float(line[5])/float(line[2])) > 0.9 and line[1] != line[3] and line[2] != line[5] and float(float(line[4])) > 93.0:
        #quality_hits.append(line)
#        print("allele", line)
        allele_hits.append(line)
    elif float(float(line[5])/float(line[2])) > 0.9 and line[1] != line[3] and line[2] != line[5] and float(float(line[4])) < 93.0 and float(float(line[4])) > 75.0:
        ohno_hits.append(line)
#        print("ohnologue", line)
    elif float(float(line[5])/float(line[2])) == 1.0 and line[1] == line[3] and line[2] == line[5]:
        hits_to_self.append(line)

allele_names = []
for line in allele_hits:
    allele_names.append(line[0])

gene_dict = {}
for gene in genes:
    g = gene.strip().split("\t")
    scaff = g[0].strip("av")
    new_scaff = 'scaffold_'+scaff
    gene_dict[g[1]] = [new_scaff, g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]]



for gene in gene_dict.keys():
    print gene_dict[gene]

def check_for_gene_overlap(unpaired_list, dict): # This function will search through gene_dict for a matching gene regrion to the quality_hits line
    for key in dict.keys():
        if dict[key][0] == unpaired_list[3]: # Look for the matching scaffold to the line
            a = range(int(unpaired_list[10]), int(unpaired_list[11]))
            b = range(int(dict[key][3]), int(dict[key][4]))
            a_s = set(a)
            new_gene = [] # this is a list holder that we can store relevant information if we find a match to an existing gene.
            if len(a_s.intersection(b)) > 1: # Im starting with one here, but I don't kow what constitutes an already existing gene. Probably more than the ovelap of one pos#ition.
                
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
                new_gene.append([unpaired_list[0], unpaired_list[1],new_name,dict[key][0],unpaired_list[10],unpaired_list[11],dict[key][3],dict[key][4], unpaired_list[4]])
#                print "we found gene %s scaffold %s with a start %d and stop %d that matches %s scaffold with a start %s and end %s" %(gene_dict[key][2],gene_dict[key][0],int(gene_dict[key][3]),int(gene_dict[key][4]),unpaired_list[3], unpaired_list[10], unpaired_list[11])
                return new_gene

for p in allele_hits:
    w = check_for_gene_overlap(p, gene_dict)
    


new_genes = []
existing_genes = []
for p in allele_hits:
    w = check_for_gene_overlap(p, gene_dict)
    if w != None:
        if w[0] in existing_genes:
#            print("already in the list - a duplicate hit")
        else:
            existing_genes.append(w[0])

        
    if w == None: # These are the genes that have no overlap with existing genes.  Some of the genes within the quality_hits occur more than once, because they have multiple quality hits to other scaffolds.  This script record show many other hits there are
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
        new_genes.append([p[0],p[1], new_name, p[3], p[10], p[11]])
