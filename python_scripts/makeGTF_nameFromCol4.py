# Modified 10/19/21 to make the GTF 1-based instead of 0-based

import sys

inBed = open(sys.argv[1],'rt')
outGTF = open(sys.argv[2],'wt')

count=1
for line in inBed:
	splitLine = line.strip().split("\t")
	outGTF.write(str(splitLine[0])+"\t"+"peakCaller"+"\t"+"exon"+"\t"+str(int(splitLine[1])+1)+"\t"+str(splitLine[2])+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"gene_id \""+str(splitLine[3])+"\""+";"+"\n")
	count+=1
inBed.close()
outGTF.close()
