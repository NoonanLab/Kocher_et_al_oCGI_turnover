# 7/26/22
# Purpose: take bed file and sort names in each field in column 4
# Where names are a_###### or b_######, or a mix of the two types
# Sort as such:
# 'a' first (region # ascending), then 'b' (region # ascending)

# python sortRegionNames.py <bed>

import sys

inBed = open(sys.argv[1],'rt')

for line in inBed:
    splitLine = line.strip().split()
    regionNames = splitLine[3].split(',')

    aRegions = []
    bRegions = []

    finalRegionNames = []
    
    for i in regionNames:
        if i[0] == 'a':
            aRegions.append(int(i.split('_')[1]))
        if i[0] == 'b':
            bRegions.append(int(i.split('_')[1]))

    for i in sorted(aRegions):
        outString = 'a_'+str(i)
        finalRegionNames.append(outString)
    for i in sorted(bRegions):
        outString = 'b_'+str(i)
        finalRegionNames.append(outString)

    print(str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[2])+'\t'+','.join(finalRegionNames))

inBed.close()
        
