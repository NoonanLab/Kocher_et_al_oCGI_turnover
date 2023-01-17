# 7/30/22
# Purpose: take count files (output by featureCounts) for each replicate and output the average RPM, only for reconciledPeak regions
# Updated 9/22/22 to remove step filtering for regions starting with a_ or b_ to work in singleSpecies pipeline instead of speciesPairs pipeline

import sys

inFileNumber = int(sys.argv[1])

fileArray = []

TotalReads = []
# TotalReads = [# of reads in file 1, # of reads in file 2...]

if inFileNumber == 1:
    fileArray.append(sys.argv[2])
    TotalReads = [0]
if inFileNumber == 2:
    fileArray.append(sys.argv[2])
    fileArray.append(sys.argv[3])
    TotalReads = [0,0]
elif inFileNumber == 3:
    fileArray.append(sys.argv[2])
    fileArray.append(sys.argv[3])
    fileArray.append(sys.argv[4])
    TotalReads = [0,0,0]
elif inFileNumber == 4:
    fileArray.append(sys.argv[2])
    fileArray.append(sys.argv[3])
    fileArray.append(sys.argv[4])
    fileArray.append(sys.argv[5])
    TotalReads = [0,0,0,0]

Counts_Dict = {}
# Counts_Dict[regionName] = [counts in rep 1, counts in rep 2...]

Length_Dict = {}
# Length_Dict[regionName] = # bp

index = 0
for fileName in fileArray:
    
    inRepFile = open(fileName,'rt')
    for line in inRepFile:
        if line[0] != '#' and line[0] != 'G':
            splitLine = line.strip().split()
            regionName = splitLine[0]
            chrom = splitLine[1]
            start = splitLine[2]
            end = splitLine[3]
            length = float(splitLine[5]) # this is corrected for the gtf being 1-based
            readCount = float(splitLine[6])

            TotalReads[index] += readCount
            if regionName not in Counts_Dict:
                Counts_Dict[regionName] = ['.'] * inFileNumber
                Length_Dict[regionName] = length
            Counts_Dict[regionName][index] = readCount
    inRepFile.close()
    index += 1

# TotalReads = [# of reads in file 1, # of reads in file 2...]
# Counts_Dict[regionName] = [counts in rep 1, counts in rep 2...]
# Length_Dict[regionName] = # bp
    
for regionName in Counts_Dict:
    #print(regionName)
    #print(Counts_Dict[regionName])
    #print(TotalReads)
    RPMarray = []
    RPKMarray = []
    for i in range(0,inFileNumber):
        RPMarray.append(Counts_Dict[regionName][i] / TotalReads[i] * 1000000)
        RPKMarray.append(Counts_Dict[regionName][i] / TotalReads[i] / Length_Dict[regionName] * 1000000 * 1000)
    avgRPM = sum(RPMarray) / len(RPMarray)
    avgRPKM = sum(RPKMarray) / len(RPKMarray)

    #print(avgRPM)

    print(regionName+'\t'+str(round(avgRPM,2))+'\t'+str(round(avgRPKM,2)))
    
