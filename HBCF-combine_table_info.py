#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='''Combine the information contained in the genes.all table found in JOES_CONSORTIUM with the pairs.alleles in the Avaga Consortium file and the ALLELES_NEW_GENES and ALLELES_EXISTING_GENES tables.  All of this information eventually needs to be put into a db but I need info today so I'm wirting this script''')
parser.add_argument('-genes_all', help = 'the genes.all file in JOES_CONSORTIUM which contains all genes from the Avaga CONSORTIUM and those identified during gene and pair mining', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/genes.all')
parser.add_argument('-av_pairs_alleles', help = 'the pairs.alleles file found in the Avaga consortium directory', default = '/groups/rotifer/Avgenome/Genoscope/v2/ConsortiumFiles/pairs.alleles')
parser.add_argument('-new_alleles', help = 'the new alleles identified during gene and pair mining', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/L75_A90_O70/ALLELES_NEW_GENES.txt')
parser.add_argument('-exist_alleles', help = 'the existing alleles identified during the gene and pair mining', default = '/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/GENE_AND_PAIR_MINING/L75_A90_O70/ALLELES_EXISTING_GENES.txt')
args = parser.parse_args()

# create a function for dictionary creation

def dict_creator(infile):
    d = {}
    dat = open(infile, 'rU')
    count = 0
    for line in dat:
        x = line.strip().split('\t')
        d[count] = x[0:len(x)]
        count += 1
    return d

av_pairs = dict_creator(args.av_pairs_alleles)
new_pairs = dict_creator(args.new_alleles)
exist_pairs = dict_creator(args.exist_alleles)
genes = dict_creator(args.genes_all)

outfile = open('/workspace/markwelchlab/Haplotype_Based_Conversion_Finder/JOES_CONSORTIUM_FILES/gene_pairs_with_scaffold_info.txt', 'w')

for key in av_pairs.keys():
    p1 = av_pairs[key][2]
    p2 = av_pairs[key][3]
    for key in genes.keys():
        if p1 == genes[key][1]:
            outfile.write(p1+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\t')
            print("just wrote1", p1, genes[key][0], genes[key][2], genes[key][3])
    for key in genes.keys():
        if p2 == genes[key][1]:
            outfile.write(p2+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\n')
            print("just wrote2", p2, genes[key][0], genes[key][2], genes[key][3])
for key in exist_pairs.keys():
    p1 = exist_pairs[key][0]
    p2 = exist_pairs[key][2]
    p2_test = exist_pairs[key][2]+'_'+exist_pairs[key][4]+'_'+exist_pairs[key][5]
#    print p2_test
    for key in genes.keys():
        if p1 == genes[key][1]:
            outfile.write(p1+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\t')
    for key in genes.keys():
        if p2_test == genes[key][1]+'_'+genes[key][2]+'_'+genes[key][3]:
#            print p2, genes[key][0], genes[key][2], genes[key][3]
            outfile.write(p2+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\n')

for key in new_pairs.keys():
    p1 = new_pairs[key][0]
    p2 = new_pairs[key][2]
    for key in genes.keys():
        if p1 == genes[key][1]:
            outfile.write(p1+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\t')
    for key in genes.keys():
        if p2 == genes[key][1]:
            outfile.write(p2+'\t'+genes[key][0]+'\t'+str(genes[key][2])+'\t'+str(genes[key][3])+'\n')

outfile.close()
