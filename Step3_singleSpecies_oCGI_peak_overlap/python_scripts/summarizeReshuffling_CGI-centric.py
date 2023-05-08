# 1/3/23
# Purpose: take results of reshuffling and output obs value, exp value, and p value

import sys
import glob

inSpecies = sys.argv[1]
inTissue = sys.argv[2]
inMark = sys.argv[3]

inSpeciesCounts = inSpecies
inMarkCounts = inMark

# rename a second species & mark variable to match formatting in countsWithPeak/ files
if inTissue == 'devBrain' or inTissue == 'devLimb':
    if inSpecies == 'rheMac2':
        inSpeciesCounts = 'rheMac10'
    if inSpecies == 'mm9':
        inSpeciesCounts = 'mm39'
        
    if inMark == 'ac':
        inMarkCounts = 'H3K27ac'
    if inMark == 'me2':
        inMarkCounts = 'H3K4me2'

inSpeciesSummary = open('observedCounts.txt', 'rt')
inCountFileList = glob.glob('countsWithPeak/'+inSpeciesCounts+'_'+inTissue+'_'+inMarkCounts+'_countsWithPeak_jobNum*.txt')

# get counts for total oCGIs and observed oCGIs with a peak for this tissue and mark
for line in inSpeciesSummary:
    splitLine = line.strip().split('\t')
    if splitLine[0] == inSpecies and splitLine[1] == inTissue and splitLine[2] == inMark:
        observedWithPeak = float(splitLine[3])
        observedDenominator = float(splitLine[4])
inSpeciesSummary.close()

# calculate observed ratio
observedRatio = observedWithPeak / observedDenominator

# initialize variables for storing info on reshuffling
countWhereShuffledIsHigher = 0
countWhereShuffledIsLower = 0
ratioArray = []

# read results of reshuffling test
for file in inCountFileList:
    inCountFile = open(file, 'rt')
    
    for line in inCountFile:
    
        # compare observed to reshuffled values
        splitLine = line.strip().split('\t')
        shuffledWithPeak = float(splitLine[0])
        shuffledDenominator = float(splitLine[1])
        
        # calculate shuffled ratio (because not all CGIs are able to reshuffle)
        shuffledRatio = shuffledWithPeak / shuffledDenominator
    
        # see if higher or lower than observedRatio
        if shuffledRatio > observedRatio:
           countWhereShuffledIsHigher += 1
        if shuffledRatio < observedRatio:
            countWhereShuffledIsLower += 1
        
        # append ratio to ratioArray for getting expected value
        ratioArray.append(shuffledRatio)
        
    inCountFile.close()
    
# calculate expected value using average of ratioArray
expectedRatio = sum(ratioArray) / len(ratioArray)

# calculate two-sided p-value
#print(countWhereShuffledIsHigher)
#print(countWhereShuffledIsLower)

if countWhereShuffledIsHigher < countWhereShuffledIsLower:
    if countWhereShuffledIsHigher == 0:
        countWhereShuffledIsHigher = 1
    p = float(countWhereShuffledIsHigher) / 20000
if countWhereShuffledIsHigher > countWhereShuffledIsLower:
    if countWhereShuffledIsLower == 0:
        countWhereShuffledIsLower = 1
    p = float(countWhereShuffledIsLower) / 20000
     
# output results
print(inSpeciesCounts +'\t'+ inTissue +'\t'+ inMarkCounts +'\t'+ str(observedRatio) +'\t'+ str(expectedRatio) +'\t' + str(p))    

