# 2/23/22
# Purpose: count number of peaks for each sample in the following categories:
# 1) peaks that start before hs754 and end within it
# 2) peaks that start before hs754 and end downstream of it by > 390bp
# 2.5) peaks that start before hs54 and end downstream but by less than 390bp)
# 3) peaks contained fully within hs754 (i.e. 5531bp after the end of the mouse 5' homology arm)
# 4) peaks that start within hs754 and end downstream of it
# 5) peaks that start within 390bp downstream of hs754 - should they be shifted back the full 390bp?
# 6) total peaks on chr13

# Usage: python summarizePeakCategories.py <*.narrowPeak/broadPeak>

import sys

inPeaks = open(sys.argv[1],'rt')
sampleString = '_'.join(sys.argv[1].split('_')[0:6])

# store numbers of peaks in each category as you go through the lines
category1 = 0
category2 = 0
category2_5 = 0
category3 = 0
category4 = 0
category5 = 0
category6 = 0

for line in inPeaks:
    splitLine = line.strip().split()
    chrom = splitLine[0]
    start = int(splitLine[1])
    end = int(splitLine[2])

    if chrom == 'chr13':
        # add to total peaks
        category6 += 1

        # add to all other categories based on locations of start and end
        if start < 72436073:
            if end > 72436073 and end <= (72436073 + 5141):
                category1 += 1
            elif end > (72436073 + 5141 + 390):
                category2 += 1
            elif end > (72436073 + 5141) and end < (72436073 + 5141 + 390):
                category2_5 += 1
        elif start >= 72436073 and start < (72436073 + 5141):
            if end <= (72436073 + 5141):
                category3 += 1
            elif end > (72436073 + 5141):
                category4 += 1
        elif start >= (72436073 + 5141) and start <= (72436073 + 5141 + 390):
            category5 += 1
inPeaks.close()
print(sampleString+'\t'+str(category1)+'\t'+str(category2)+'\t'+str(category2_5)+'\t'+str(category3)+'\t'+str(category4)+'\t'+str(category5)+'\t'+str(category6))
    

    
