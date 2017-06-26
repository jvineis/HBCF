#!/bin/usr/env python

import sys 

infile = open(sys.argv[1], 'rU')
outfile = open(sys.argv[2], 'w')
contig_id = sys.argv[3]
outfile.write("gene_callers_id"+"\t"+"contig"+"\t"+"start"+"\t"+"stop"+"\t"+"direction"+"\t"+"partial"+"\t"+"source"+"\t"+"version"+"\n")

num = 0
for line in infile:
    x = line.strip().split("\t")
    if int(x[2]) > int(x[3]):
        outfile.write(str(num)+"\t"+contig_id+"\t"+x[3]+"\t"+x[2]+"\t"+"f"+"\t"+"1"+"\t"+"Consortium"+"\t"+"v2"+"\n")
    elif int(x[2]) < int(x[3]):
        outfile.write(str(num)+"\t"+contig_id+"\t"+x[2]+"\t"+x[3]+"\t"+"f"+"\t"+"1"+"\t"+"Consortium"+"\t"+"v2"+"\n")
    num += 1
    
    
