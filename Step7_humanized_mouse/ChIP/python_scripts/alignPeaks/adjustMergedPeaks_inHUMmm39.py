# 2/23/22
# Purpose: shifts peaks in MERGED peak files to match where they should be on HUM mm39

# Usage: python adjustHUMpeaksInMm39.py <merged peak file in mm39 coordinates> <${sampleName}_mergedPeaks_renamed_fallingInHs754_LOhg38_backToHUM.bed (peaks falling totally in hs754 lifted from mm39 coordinates to hg38, then into HUMmm39)> > <${sampleName}_mergedPeaks_inHUMmm39.bed (output file for use in HTSeq)>
# python adjustMergedPeaks_inHUMmm39.py ${i}_mergedPeaks_renamed.bed ${sampleName}_mergedPeaks_renamed_fallingInHs754_LOhg38_backToHUM.bed > ${sampleName}_mergedPeaks_inHUMmm39.bed

import sys

inMergedPeaksInMm39 = open(sys.argv[1],'rt')
inHs754peaksInHUMmm39 = open(sys.argv[2],'rt')

# store peaks falling completely within hs754, moved to hg38 coord and then lifted to mouse
hs754_List = []
for line in inHs754peaksInHUMmm39:
    hs754_List.append(line.strip().split(';')[0])
inHs754peaksInHUMmm39.close()

haveHs754PeaksBeenPrinted = 'no'

# go through file of peaks in HUM coordinates and print to output the entries adjusted into mm39
for line in inMergedPeaksInMm39:
    splitLine = line.strip().split()
    chrom = splitLine[0]
    start = int(splitLine[1])
    end = int(splitLine[2])
    name = splitLine[3].split(';')[0]

    if chrom == 'chr13':
        if start < 72436073 and end < 72436073:
            print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)
        elif start < 72436073 and end > 72441604:
            newEnd = end + 390
            print(chrom+'\t'+str(start)+'\t'+str(newEnd)+'\t'+name)
        elif start > 72436073 and start < 72441604 and end > 72436073 and end < 72441604:
            if haveHs754PeaksBeenPrinted == 'no':
                print('\n'.join(hs754_List))
                haveHs754PeaksBeenPrinted = 'yes'
        # peaks that start within hs754 and end downstream
        # only the case for D_e17.5_H3K27ac_HUM which was already printed above - so do nothing
        elif start > 72436073 and start < 72441604 and end > 72436073:
            pass
        elif start > 72441604 and end > 72441604:
            newStart = start + 390
            newEnd = end + 390
            print(chrom+'\t'+str(newStart)+'\t'+str(newEnd)+'\t'+name)
        else:
            print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)

    else:
        print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+name)
inMergedPeaksInMm39.close()
