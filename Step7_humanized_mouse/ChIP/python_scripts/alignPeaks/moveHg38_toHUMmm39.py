# 2/23/22
# Purpose: take peaks in hg38 coordinates and move them to HUM mm39 coordinates
# have to also reverse direction

import sys

inPeaks = open(sys.argv[1],'rt')

mouseRegionChrom = 'chr13'
humanRegionStart = 3193946
humanRegionEnd = 3199477

mouseRegionStart = 72436073
mouseRegionEnd = 72441214

for line in inPeaks:
    splitLine = line.strip().split()
    peakInHumanStart = int(splitLine[1])
    peakInHumanEnd = int(splitLine[2])

    startInHUMmouseCoord = mouseRegionStart + (humanRegionEnd - peakInHumanEnd)
    endInHUMmouseCoord = mouseRegionEnd + 390 - (peakInHumanStart - humanRegionStart)

    print(mouseRegionChrom+'\t'+str(startInHUMmouseCoord)+'\t'+str(endInHUMmouseCoord)+'\t'+str(splitLine[3]))    #+'\t'+'\t'.join(splitLine[3:]))

inPeaks.close()
