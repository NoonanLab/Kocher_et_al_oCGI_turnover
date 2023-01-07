# 7/31/22
# Pipeline in which this is used:
# Start with bed file in species X (bed #1), with unique names for each row in column 4
# Lift bed #1 to species Y (bed #2)
# Then lift bed #2 back to species X (bed #3)
# Use bedtools intersect -wao to intersect bed #3 with bed #1
# This script then takes that intersection output and writes ONLY those sites that intersect themselves, in species X coordinates (bed #4)

import sys

inIntersection = open(sys.argv[1],'rt')

Overlap_Dict = {}
# Overlap_Dict[liftedBackName] = [intersectInOriginalName 1, intersectInOriginalName 2, ...]

Coord_Dict = {}
# Coord_Dict[liftedBackName] = [chr,start,end]

UniqueList = []

for line in inIntersection:
    splitLine = line.strip().split()
    chrom = splitLine[0]
    start = splitLine[1]
    end = splitLine[2]
    liftedBackName = splitLine[3]
    intersectInOriginalName = splitLine[7]
    
    if liftedBackName not in Overlap_Dict:
        Overlap_Dict[liftedBackName] = []
        Coord_Dict[liftedBackName] = [chrom,start,end]

    if intersectInOriginalName != '.':
        Overlap_Dict[liftedBackName].append(intersectInOriginalName)
inIntersection.close()

for liftedBackName in Overlap_Dict:
    if len(Overlap_Dict[liftedBackName]) == 1:
        print('\t'.join(Coord_Dict[liftedBackName])+'\t'+Overlap_Dict[liftedBackName][0])
