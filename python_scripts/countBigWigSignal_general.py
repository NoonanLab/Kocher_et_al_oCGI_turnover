# 8/10/22
# Purpose: take count files (output by bigWigAverageOverBed) for each replicate and output the average counts, only for reconciledPeak regions

import sys

inFileNumber = int(sys.argv[1])

#The output columns are:
#   name - name field from bed, which should be unique
#   size - size of bed (sum of exon sizes
#   covered - # bases within exons covered by bigWig
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases

fileArray = []

TotalSignal = []
# TotalSignal = [sum of signal in file 1, sum of signal in file 2...]

if inFileNumber == 1:
    fileArray.append(sys.argv[2])
    TotalSignal = [0]
if inFileNumber == 2:
    fileArray.append(sys.argv[2])
    fileArray.append(sys.argv[3])
    TotalSignal = [0,0]
elif inFileNumber == 3:
    fileArray.append(sys.argv[2])
    fileArray.append(sys.argv[3])
    fileArray.append(sys.argv[4])
    TotalSignal = [0,0,0]

Signal_Dict = {}
# Signal_Dict[regionName] = [total signal in rep 1, total signal in rep 2...]

Length_Dict = {}
# Length_Dict[regionName] = # bp

index = 0
for fileName in fileArray:
    
    inRepFile = open(fileName,'rt')
    for line in inRepFile:
        splitLine = line.strip().split()
        regionName = splitLine[0]
        length = float(splitLine[1])
        signal = float(splitLine[3])

        TotalSignal[index] += signal
        if regionName not in Signal_Dict:
            Signal_Dict[regionName] = ['.'] * inFileNumber
            Length_Dict[regionName] = length
        Signal_Dict[regionName][index] = signal
    inRepFile.close()
    index += 1

# TotalSignal = [sum of signal in file 1, sum of signal in file 2...]
# Signal_Dict[regionName] = [total signal in rep 1, total signal in rep 2...]
# Length_Dict[regionName] = # bp
    
for regionName in Signal_Dict:
    #print(regionName)
    #print(Counts_Dict[regionName])
    #print(TotalReads)
    RPMarray = []
    RPKMarray = []
    for i in range(0,inFileNumber):
        RPMarray.append(Signal_Dict[regionName][i] / TotalSignal[i] * 1000000)
        RPKMarray.append(Signal_Dict[regionName][i] / TotalSignal[i] / Length_Dict[regionName] * 1000000 * 1000)
    avgRPM = sum(RPMarray) / len(RPMarray)
    avgRPKM = sum(RPKMarray) / len(RPKMarray)

    #print(avgRPM)

    print(regionName+'\t'+str(round(avgRPM,2))+'\t'+str(round(avgRPKM,2)))
    
