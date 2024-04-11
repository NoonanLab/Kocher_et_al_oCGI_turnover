# 3/5/3024
# Purpose: take 3 bed files and make output bed files that only contain lines for chromosomes present in all 3 files

import sys

inFile1 = open(sys.argv[1], 'rt')
inFile2 = open(sys.argv[2], 'rt')
inFile3 = open(sys.argv[3], 'rt')

outFile1 = open(sys.argv[4], 'wt')
outFile2 = open(sys.argv[5], 'wt')
outFile3 = open(sys.argv[6], 'wt')

File1_List = []
File2_List = []
File3_List = []

for line in inFile1:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom not in File1_List:
        File1_List.append(chrom)
inFile1.close()

for line in inFile2:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom not in File2_List:
        File2_List.append(chrom)
inFile2.close()

for line in inFile3:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom not in File3_List:
        File3_List.append(chrom)
inFile3.close()

#print(len(File1_List))
#print(File1_List)
#print(File2_List)
#print(File3_List)

# get chromosomes present in all 3 lists

commonChrom_List = []

for chrom in File1_List:
    #print(chrom)
    if chrom in File2_List:
        if chrom in File3_List:
            commonChrom_List.append(chrom)
        
#print(commonChrom_List)
        
# print to output files

inFile1 = open(sys.argv[1], 'rt')
inFile2 = open(sys.argv[2], 'rt')
inFile3 = open(sys.argv[3], 'rt')

for line in inFile1:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom in commonChrom_List:
        outFile1.write(line)
inFile1.close()
outFile1.close()

for line in inFile2:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom in commonChrom_List:
        outFile2.write(line)
inFile2.close()
outFile2.close()

for line in inFile3:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    if chrom in commonChrom_List:
        outFile3.write(line)
inFile3.close()
outFile3.close()
