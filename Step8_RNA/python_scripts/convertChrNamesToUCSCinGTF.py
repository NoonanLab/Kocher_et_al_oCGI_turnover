# 11/21/22
# Purpose: take ENSEMBL GTF file and convert chromosome names to UCSC format
# Step 1: test for chromosome number of 1-2 characters
# Step 2: add 'chr' before number and print to output

import sys

inGTF = open(sys.argv[1], 'rt')
outGTF = open(sys.argv[2], 'wt')

for line in inGTF:
    if line[0] == '#':
        outGTF.write(line.strip() + '\n')
    else:
        splitLine = line.strip().split('\t')
        chrom = str(splitLine[0])
        if len(chrom) < 3:
            newChrom = 'chr'+chrom
            splitLine[0] = newChrom
            outGTF.write('\t'.join(splitLine) + '\n')
inGTF.close()
outGTF.close()