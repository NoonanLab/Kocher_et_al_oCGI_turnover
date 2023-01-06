# 7/27/22
# Purpose: take peak numbers and output bedtools commands for making intersection files based on peak numbers

# python makeIntersectCommands.py peakSummary.txt repsToUseRoller.txt > intersectCommands.txt

import sys

inPeakSummary = open(sys.argv[1],'rt')
outRepsToUse = open(sys.argv[2],'wt')

outputArray = []

for line in inPeakSummary:
    splitLine = line.strip().split()
    species = splitLine[0]
    tissue = splitLine[1]
    mark = splitLine[2]
    rep1 = float(splitLine[3])
    rep2 = float(splitLine[4])
    rep3 = float(splitLine[5])

    # get average of each combination of 2 replicates, and store in an array
    avgArray = [(rep1+rep2)/2, (rep1+rep3)/2, (rep2+rep3)/2]

    # store max average
    maxAvg = max(avgArray)

    # store difference of each replicate to the max number of peaks
    diffArray = [(rep1-maxAvg)/maxAvg*100, (rep2-maxAvg)/maxAvg*100, (rep3-maxAvg)/maxAvg*100]

    repToUseArray = [1,2,3]
    
    # remove a replicate if it has 50%+ more peaks or 20%+ fewer peaks than average of other two replicates
    for rep in [1,2,3]:
        if diffArray[rep-1] >= 50:
            repToUseArray.remove(rep)

    if len(repToUseArray) == 3:
        for rep in [1,2,3]:
            if diffArray[rep-1] <= -20:
                repToUseArray.remove(rep)
    
    # generate intersection commands
    if mark == 'H3K4me3' or mark == 'H3K27ac':
        peakType = 'narrow'
    elif mark == 'H3K4me1':
        peakType = 'broad'
    
    if len(repToUseArray) == 3:
        outputArray.append('bedtools intersect -a '+species+'_'+tissue+'_'+mark+'_1_peaks.'+peakType+'Peak -b '+species+'_'+tissue+'_'+mark+'_2_peaks.'+peakType+'Peak | bedtools intersect -a - -b '+species+'_'+tissue+'_'+mark+'_3_peaks.'+peakType+'Peak | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\t'+species+'_'+tissue+'_'+mark+'_\"NR }\' > intersection/'+species+'_'+tissue+'_'+mark+'_intersection.bed')
    elif len(repToUseArray) == 2:
        outputArray.append('bedtools intersect -a '+species+'_'+tissue+'_'+mark+'_'+str(repToUseArray[0])+'_peaks.'+peakType+'Peak -b '+species+'_'+tissue+'_'+mark+'_'+str(repToUseArray[1])+'_peaks.'+peakType+'Peak | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\t'+species+'_'+tissue+'_'+mark+'_\"NR }\' > intersection/'+species+'_'+tissue+'_'+mark+'_intersection.bed')

    stringArray = []
    for i in repToUseArray:
        stringArray.append(str(i))
    outRepsToUse.write(species+'\t'+tissue+'\t'+mark+'\t'+','.join(stringArray)+'\n')

inPeakSummary.close()
outRepsToUse.close()

# print intersection commands
print('cd /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; '+' ; '.join(outputArray))
