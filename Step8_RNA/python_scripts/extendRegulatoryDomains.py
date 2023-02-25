# 12/12/22
# Purpose: take bed file with promoters (5kb upstream and 1kb downstream of TSSs)
# and extend them to "regulatory domains" = extend to nearest one or up to 1 Mb

import sys

inPromoters = open(sys.argv[1], 'rt')
inChromSizes = open(sys.argv[2], 'rt')

ChromSizes_Dict = {}
# ChromSizes_Dict[chr] = length

for line in inChromSizes:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    length = int(splitLine[1])
    ChromSizes_Dict[chrom] = length

Coord_Dict = {}
# Coord_Dict[chrom][n] = [start, end, ENSG]

previousChrom = '.'
chromList = []

for line in inPromoters:

    # get current coordinates
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    start = int(splitLine[1])
    end = int(splitLine[2])
    ENSG = splitLine[3]
    
    # start new counter and initialize new sub-dictionary if necessary
    if chrom != previousChrom:
        n = 1
        Coord_Dict[chrom] = {}
        chromList.append(chrom)

    # store in Coord_Dict
    Coord_Dict[chrom][n] = [start, end, ENSG]
        
    n += 1
    previousChrom = chrom
inPromoters.close()

#for n in Coord_Dict['chr1']:
#    print(str(n))
#    print(Coord_Dict['chr1'][n])

# loop through all chromosomes
for chrom in chromList:
    # write to output for first entry
    n = 1
    start = Coord_Dict[chrom][n][0]
    end = Coord_Dict[chrom][n][1]
    ENSG = Coord_Dict[chrom][n][2]
    
    nextStart = Coord_Dict[chrom][n + 1][0]
    
    distanceToZero = start - 0
    distanceToNextStart = nextStart - start
    
    # get new start
    if distanceToZero <= 1000000:
        newStart = 0
    elif distanceToZero > 1000000:
        newStart = start - 1000000
        
    # get new end
    if nextStart <= end:
        newEnd = end
    
    elif nextStart > end:
        if distanceToNextStart <= 1000000:
            newEnd = nextStart
        elif distanceToNextStart > 1000000:
            newEnd = end + 1000000
    
    print(chrom +'\t'+ str(newStart) +'\t'+ str(newEnd) +'\t'+ ENSG)
    
    # write to output for all middle entries
    n = 2
    while n < max(Coord_Dict[chrom].keys()):
    
        start = Coord_Dict[chrom][n][0]
        end = Coord_Dict[chrom][n][1]
        ENSG = Coord_Dict[chrom][n][2]
        previousEnd = Coord_Dict[chrom][n - 1][1]
        nextStart = Coord_Dict[chrom][n + 1][0]
        distanceToPreviousEnd = start - previousEnd
        distanceToNextStart = nextStart - start
    
        # get new start
        if start <= previousEnd:
            newStart = start
        
        elif start > previousEnd:
            if distanceToPreviousEnd <= 1000000:
                newStart = previousEnd
            elif distanceToPreviousEnd > 1000000:
                newStart = start - 1000000
    
        # get new end
        if end >= nextStart:
            newEnd = end
            
        elif end < nextStart:
            if distanceToNextStart <= 1000000:
                newEnd = nextStart
            elif distanceToNextStart > 1000000:
                newEnd = end + 1000000
        
        # print
        print(chrom +'\t'+ str(newStart) +'\t'+ str(newEnd) +'\t'+ ENSG)
        
        # increase n
        n += 1
    
    # write to output for last entry
    n = max(Coord_Dict[chrom].keys())
    start = Coord_Dict[chrom][n][0]
    end = Coord_Dict[chrom][n][1]
    ENSG = Coord_Dict[chrom][n][2]

    # get new start
    if start <= previousEnd:
        newStart = start
        
    elif start > previousEnd:    
        if distanceToPreviousEnd <= 1000000:
            newStart = previousEnd
        elif distanceToPreviousEnd > 1000000:
            newStart = start - 1000000
    
    # get new end
    chromEnd = ChromSizes_Dict[chrom]
    distanceToChromEnd = chromEnd - end

    if distanceToChromEnd <= 1000000:
        newEnd = chromEnd
    elif distanceToChromEnd > 1000000:
        newEnd = end + 1000000
    
    # print
    print(chrom +'\t'+ str(newStart) +'\t'+ str(newEnd) +'\t'+ ENSG)
