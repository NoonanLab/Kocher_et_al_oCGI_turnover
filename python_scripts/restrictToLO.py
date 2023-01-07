# 3/29/21
# Purpose: reprint lines in a bed file (file A), but only if col 4 is contained in another bed file (file B)
# to make a new output (file C)

# usage: python restrictToLO.py fileA.bed fileB.bed > fileC.bed

import sys

inFileA = open(sys.argv[1],'rt')
inFileB = open(sys.argv[2],'rt')

peakLabelList = []
for line in inFileB:
    splitLine = line.strip().split()
    if splitLine[3] not in peakLabelList:
        peakLabelList.append(splitLine[3])
inFileB.close()

for line in inFileA:
    splitLine = line.strip().split()
    if splitLine[3] in peakLabelList:
        print(str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[2])+'\t'+str(splitLine[3]))
inFileA.close()
