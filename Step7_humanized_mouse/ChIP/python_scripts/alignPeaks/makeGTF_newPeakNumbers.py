# 2/23/22
# Purpose: turn bed file into GTF, and in the process output a new bedfile with peak numbers in addition to names. Use new peak numbers in GTF.

import sys

inBed = open(sys.argv[1],'rt')
outGTF = open(sys.argv[2],'wt')
outBed = open(sys.argv[3],'wt')

#count=1
for line in inBed:
    splitLine = line.strip().split()
    name = splitLine[0]+':'+splitLine[1]+'-'+splitLine[2]
    outGTF.write(str(splitLine[0])+"\t"+"peakCaller"+"\t"+"exon"+"\t"+str(splitLine[1])+"\t"+str(splitLine[2])+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+"gene_id \""+name+"\""+";"+"\n")
    outBed.write(str(splitLine[0])+"\t"+str(splitLine[1])+"\t"+str(splitLine[2])+"\t"+name+"\n")
    #count+=1
inBed.close()
outGTF.close()
outBed.close()
