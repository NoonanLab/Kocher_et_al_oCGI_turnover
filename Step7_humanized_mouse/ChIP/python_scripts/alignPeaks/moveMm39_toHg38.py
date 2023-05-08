# 2/23/22
# Purpose: take peaks in mm39_humanized coordinates and move them to hg38 coordinates
# have to also reverse direction

import sys

inPeaks = open(sys.argv[1],'rt')

humanRegionChrom = 'chr5'
humanRegionStart = 3193946
humanRegionEnd = 3199477

mouseRegionStart = 72436073
mouseRegionEnd = 72441214

for line in inPeaks:
    splitLine = line.strip().split()
    humanizedMousePeakStart = int(splitLine[1])
    humanizedMousePeakEnd = int(splitLine[2])

    startInHumanCoord = humanRegionStart + (mouseRegionEnd + 390 - humanizedMousePeakEnd)
    endInHumanCoord = humanRegionEnd - (humanizedMousePeakStart - mouseRegionStart)

    print(humanRegionChrom+'\t'+str(startInHumanCoord)+'\t'+str(endInHumanCoord)+'\t'+str(splitLine[3]))    #+'\t'+'\t'.join(splitLine[3:]))

inPeaks.close()
