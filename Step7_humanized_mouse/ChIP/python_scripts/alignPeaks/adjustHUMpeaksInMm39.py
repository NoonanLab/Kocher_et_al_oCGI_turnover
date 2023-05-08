# 2/23/22
# Purpose: shifts peaks in HUM mouse peak files to match where they should be on mm39
# modified to deal with category 4 which appears in H3K27ac data (peak starts in humanized region and ends downstream)

# Usage: python adjustHUMpeaksInMm39.py <HUM peak file in mm39_humanized coordinates> <${sampleName}_hg38LOmm39.bed (peaks falling totally in hs754 lifted from hg38 coordinates to mm39)> > <${sampleName}_HUMpeaksInMm39.bed (output file for use in merging)>

import sys

inPeaksInHumanizedMm39 = open(sys.argv[1],'rt')
inHs754peaksLiftedToMm39 = open(sys.argv[2],'rt')

# store peaks falling completely within hs754, moved to hg38 coord and then lifted to mouse
hs754_List = []
for line in inHs754peaksLiftedToMm39:
    hs754_List.append(line.strip())
inHs754peaksLiftedToMm39.close()

haveHs754PeaksBeenPrinted = 'no'

# go through file of peaks in HUM coordinates and print to output the entries adjusted into mm39
for line in inPeaksInHumanizedMm39:
    splitLine = line.strip().split()
    chrom = splitLine[0]
    start = int(splitLine[1])
    end = int(splitLine[2])
    name = splitLine[3]

    if chrom == 'chr13':
        # peaks completely upstream of hs754: just print
        if start < 72436073 and end < 72436073:
            print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)
        # peaks that start upstream and end downstream: print same start and shift end back 390bp
        elif start < 72436073 and end > 72441604:
            newEnd = end - 390
            print(chrom+'\t'+str(start)+'\t'+str(newEnd)+'\t'+name)
        # peaks that start within hs754 and end within hs754: just print coordinates that were lifted to mm39
        elif start > 72436073 and start < 72441604 and end > 72436073 and end < 72441604:
            if haveHs754PeaksBeenPrinted == 'no':
                print('\n'.join(hs754_List))
                haveHs754PeaksBeenPrinted = 'yes'
        # peaks that start within hs754 and end downstream
        # only the case for D_e17.5_H3K27ac_HUM which was already printed above - so do nothing
        elif start > 72436073 and start < 72441604  and end > 72441604:
            pass
        # peaks that start and end after hs754: shift both start and end back 390bp
        elif start > 72441604 and end > 72441604:
            newStart = start - 390
            newEnd = end - 390
            print(chrom+'\t'+str(newStart)+'\t'+str(newEnd)+'\t'+name)
        else:
            print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)

    else:
        print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)
    
    
