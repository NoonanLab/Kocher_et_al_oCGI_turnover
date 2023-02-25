# 11/22/22
# Purpose: make a bed file listing coordinates for species-specific CGIs with species-specific peaks
# for each speciesPair x tissue x mark

# INPUTS:
# 1) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles/rheMac10_mm39_brain_H3K4me3.txt 
# 2) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39/rheMac10_mm39_speciesA_reconciledCGI.bed
# 3) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39/rheMac10_mm39_speciesB_reconciledCGI.bed

# OUTPUTS:
# 1) bed file with all A-only CGIs with A-only peaks and B-only CGIs with B-only peaks, in A coordinates
# 1) bed file with all A-only CGIs with A-only peaks and B-only CGIs with B-only peaks, in B coordinates

import sys

inSummaryFile = open(sys.argv[1], 'rt')
inCGIbed_A = open(sys.argv[2], 'rt')
inCGIbed_B = open(sys.argv[3], 'rt')

speciesA = sys.argv[1].split('/')[-1].split('_')[0]
speciesB = sys.argv[1].split('/')[-1].split('_')[1]
tissue = sys.argv[1].split('/')[-1].split('_')[2]
mark = sys.argv[1].split('/')[-1].split('_')[3].split('.')[0]

outBed_A = open('nearestGene/'+speciesA+'_'+speciesB+'_'+tissue+'_'+mark+'_coordInA.bed', 'wt')
outBed_B = open('nearestGene/'+speciesA+'_'+speciesB+'_'+tissue+'_'+mark+'_coordInB.bed', 'wt')

A_CGIs_with_A_Peaks = []
B_CGIs_with_B_Peaks = []

for line in inSummaryFile:
    if line[0] != 'C':
        splitLine = line.strip().split('\t')
        
        # get CGI info
        CGIname = splitLine[0]
        CGI_A = int(splitLine[1])
        CGI_B = int(splitLine[2])
        
        # get peak info
        Peak_A = int(splitLine[18])
        Peak_B = int(splitLine[19])
        
        # if CGI & Peak are A-only, or CGI & Peak are B-only, store in appropriate list
        if CGI_A == 1 and CGI_B == 0 and Peak_A == 1 and Peak_B == 0:
            A_CGIs_with_A_Peaks.append(CGIname)
        if CGI_A == 0 and CGI_B == 1 and Peak_A == 0 and Peak_B == 1:
            B_CGIs_with_B_Peaks.append(CGIname)
            
inSummaryFile.close()

# read in coordinates from consensusCGI bed files

for line in inCGIbed_A:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    start = splitLine[1]
    end = splitLine[2]
    CGIname = splitLine[3]
    
    if CGIname in A_CGIs_with_A_Peaks or CGIname in B_CGIs_with_B_Peaks:
        outBed_A.write(chrom +'\t'+ start +'\t'+ end +'\t'+ CGIname +'\n')
inCGIbed_A.close()
    
for line in inCGIbed_B:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    start = splitLine[1]
    end = splitLine[2]
    CGIname = splitLine[3]
    
    if CGIname in A_CGIs_with_A_Peaks or CGIname in B_CGIs_with_B_Peaks:
        outBed_B.write(chrom +'\t'+ start +'\t'+ end +'\t'+ CGIname +'\n')
inCGIbed_B.close()
    