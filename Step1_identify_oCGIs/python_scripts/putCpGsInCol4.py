# 3/5/23
# Purpose: take output of faCount and output as a bed file, with number of CpGs in column 4

import sys

inFaCount = open(sys.argv[1], 'rt')
outBed = open(sys.argv[2], 'wt')

for line in inFaCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        
        # get coordinates
        coord = splitLine[0]
        
        chrom = coord.split(':')[0]
        start = coord.split(':')[1].split('-')[0]
        end = coord.split('-')[1]
        
        # get number of CpGs
        CpGs = splitLine[7]
        
        # write to output bed file
        outBed.write(chrom + '\t' + start + '\t' + end + '\t' + 'CpG: ' + CpGs + '\n')
        
inFaCount.close()
outBed.close()