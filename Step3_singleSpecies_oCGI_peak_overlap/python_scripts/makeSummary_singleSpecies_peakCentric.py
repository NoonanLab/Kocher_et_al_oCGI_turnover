# 10/14/22
# Purpose: make summary files for peakCentric singleSpecies
# IF YOU RE-RUN THIS: change header lengthPeak_Roller to lengthPeak since this script is for Noonan tissues too

import sys

# python /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/Scripts/makeSummary_singleSpecies_peakCentric.py 

inPeaks = open(sys.argv[1], 'rt')                   # rheMac10_brain_H3K27ac_noFeatures.bed 
inPeaksInHuman = open(sys.argv[2], 'rt')            # liftToHuman/rheMac10_brain_H3K27ac_peaksThatMapBack_inHumanCoord.bed 
inPeaksInHumanNoFeatures = open(sys.argv[3], 'rt')  # liftToHuman/rheMac10_brain_H3K27ac_noFeaturesInHuman_inHumanCoord.bed 
inIntersectCGIs = open(sys.argv[4], 'rt')           # intersectCGIs/rheMac10_brain_H3K27ac_intersectCGIs_RollerCoord.txt
inFeatureCounts = open(sys.argv[5], 'rt')           # featureCounts/rheMac10_brain_H3K27ac.signal
inFaCountPeaks = open(sys.argv[6], 'rt')            # faCount/rheMac10_brain_H3K27ac_Peaks.faCount
inFaCountCGIs = open(sys.argv[7], 'rt')             # faCount/rheMac10_brain_H3K27ac_CGIs.faCount
inAgeSegmentationPeaks = open(sys.argv[8], 'rt')    # intersectAge/rheMac10_brain_H3K27ac_Peaks_intersectAgeSegmentation.txt 
inAgeSegmentationCGIs = open(sys.argv[9], 'rt')     # intersectAge/rheMac10_brain_H3K27ac_CGIs_intersectAgeSegmentation.txt 
inPhastConsPeaks = open(sys.argv[10], 'rt')         # intersectPhastCons/rheMac10_brain_H3K27ac_Peaks_intersectPhastCons.txt 
inPhastConsCGIs = open(sys.argv[11], 'rt')          # intersectPhastCons/rheMac10_brain_H3K27ac_CGIs_intersectPhastCons.txt 
inRepeats = open(sys.argv[12], 'rt')                # intersectRepeats/rheMac10_brain_H3K27ac_intersectRepeats.txt

Peak_Dict = {}
# Peak_Dict[peakName]['basicInfo'] = [lengthPeak_Roller, peakLiftsToHuman, lengthPeak_Human, noFeaturesInHuman, numCpGs_Peak, GC_Peak]
# Peak_Dict[peakName]['CGI'] = [CGIname1, CGIname2, ...]
# Peak_Dict[peakName]['Peak'] = [RPM, RPKM]
# Peak_Dict[peakName]['phastCons'] = [phastCons_number, phastConsBp, phastCons_maxLOD, phastCons_sumLOD]
# Peak_Dict[peakName]['other'] = [OldestSegment_Peak, RepeatBp_Peak]

CGI_Dict = {}
# CGI_Dict[CGIname] = [lengthInRoller, NumCpGs_Roller, liftsToHuman, lengthInHuman, OldestSegment_CGI,
#                      phastConsCGI_number, phastConsCGI_bp, phastConsCGI_maxLOD, phastConsCGI_sumLOD]

AgeSegmentation_Index = {'None':0, 'Human':1, 'Ape':2, 'Primate':3, 'Eutheria':4, 'Theria':5, 'Mammalia':6, 'Amniota':7, 'Tetrapoda':8, 'Gnathostomata':9, 'Vertebrata':10}

# initialize Peak_Dict and complete Peak_Dict[peakName]['basicInfo'] = [lengthPeak_Roller, peakLiftsToHuman, lengthPeak_Human, noFeaturesInHuman, numCpGs_Peak]
for line in inPeaks:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    peakLength = int(splitLine[2]) - int(splitLine[1])
    if peakName not in Peak_Dict:
        Peak_Dict[peakName] = {}
        Peak_Dict[peakName]['basicInfo'] = [peakLength,0,0,0,0,0]
        Peak_Dict[peakName]['CGI'] = []
        Peak_Dict[peakName]['Peak'] = [0,0]
        Peak_Dict[peakName]['phastCons'] = [0,0,0,0]
        Peak_Dict[peakName]['other'] = ['None',0]
inPeaks.close()

for line in inPeaksInHuman:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    peakLengthHuman = int(splitLine[2]) - int(splitLine[1])
    Peak_Dict[peakName]['basicInfo'][1] = 1
    Peak_Dict[peakName]['basicInfo'][2] = peakLengthHuman
inPeaksInHuman.close()

for line in inPeaksInHumanNoFeatures:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    Peak_Dict[peakName]['basicInfo'][3] = 1
inPeaksInHumanNoFeatures.close()

# store info on CGI intersection in Peak_Dict[peakName]['CGI'] = [CGIname1, CGIname2, ...]
for line in inIntersectCGIs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    
    if splitLine[4] != '.':
        CGIname = splitLine[7]
        Peak_Dict[peakName]['CGI'].append(CGIname)
        lengthCGI = int(splitLine[6]) - int(splitLine[5])
        
        # also store info on CGIs in CGI_Dict[CGIname]
        if CGIname not in CGI_Dict:
            CGI_Dict[CGIname] = [lengthCGI,0,0,0,'None',0,0,0,0]
inIntersectCGIs.close()

for line in inFaCountPeaks:
    splitLine = line.strip().split('\t')
    if line[0] != '#' and line[0] != 't':
        peakName = splitLine[0].split('::')[0]
        numCpGs_Peak = int(splitLine[7])
        Peak_Dict[peakName]['basicInfo'][4] = numCpGs_Peak
        GC_content_Peak = ((float(splitLine[3]) + float(splitLine[4])) / float(splitLine[1])) * 100
        Peak_Dict[peakName]['basicInfo'][5] = GC_content_Peak
        
# store RPM / RPKM info in Peak_Dict[peakName]['Peak'] = [RPM, RPKM]
for line in inFeatureCounts:
    splitLine = line.strip().split('\t')
    peakName = splitLine[0]
    if peakName in Peak_Dict:
        RPM = splitLine[1]
        RPKM = splitLine[2]
        Peak_Dict[peakName]['Peak'][0] = RPM
        Peak_Dict[peakName]['Peak'][1] = RPKM
inFeatureCounts.close()

# store info on CGI CpGs
for line in inFaCountCGIs:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        CGIname = splitLine[0].split('::')[0]
        if CGIname in CGI_Dict:
            numCpGs_CGI = int(splitLine[7])
            CGI_Dict[CGIname][1] = numCpGs_CGI
inFaCountCGIs.close()

# store info on phastCons for Peaks in Peak_Dict[peakName]['phastCons'] = [phastCons_number, phastConsBp, phastCons_maxLOD, phastCons_sumLOD]
for line in inPhastConsPeaks:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        LOD = int(splitLine[7].split('=')[1])
        phastConsBp = int(splitLine[8])
        Peak_Dict[peakName]['phastCons'][0] += 1
        Peak_Dict[peakName]['phastCons'][1] += phastConsBp
        Peak_Dict[peakName]['phastCons'][3] += LOD
        if LOD > Peak_Dict[peakName]['phastCons'][2]:
            Peak_Dict[peakName]['phastCons'][2] = LOD

# OTHER INFO
# store info on age segmentation for Peaks
for line in inAgeSegmentationPeaks:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        if splitLine[7][0] == 'N':
            age = 'None'
        else:
            age = splitLine[7]
        overlap = int(splitLine[13])
        if overlap >= 100:
            if AgeSegmentation_Index[age] > AgeSegmentation_Index[Peak_Dict[peakName]['other'][0]]:
                Peak_Dict[peakName]['other'][0] = age
inAgeSegmentationPeaks.close()

# store info on repeat overlap for Peaks
for line in inRepeats:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    repeatOverlap = int(splitLine[8])
    Peak_Dict[peakName]['other'][1] += repeatOverlap
inRepeats.close()

# COMPLETE CGI INFO
# store info on age segmentation for CGIs
for line in inAgeSegmentationCGIs:
    splitLine = line.strip().split('\t')
    CGIname = splitLine[3]
    
    # if CGI is in this file, that means it lifts to human - also store length in human
    if CGIname in CGI_Dict:
        CGI_Dict[CGIname][2] = 1
        lengthCGI_human = int(splitLine[2]) - int(splitLine[1])
        CGI_Dict[CGIname][3] = lengthCGI_human
    
        # age segmentation
        if splitLine[4] != '.':
            if splitLine[7][0] == 'N':
                age = 'None'
            else:
                age = splitLine[7]
            overlap = int(splitLine[13])
            if overlap >= 100:
                if AgeSegmentation_Index[age] > AgeSegmentation_Index[CGI_Dict[CGIname][4]]:
                    CGI_Dict[CGIname][4] = age
inAgeSegmentationCGIs.close()

for line in inPhastConsCGIs:
    splitLine = line.strip().split('\t')
    CGIname = splitLine[3]
    if CGIname in CGI_Dict:
        if splitLine[4] != '.':
            phastConsBp = int(splitLine[8])
            LOD = int(splitLine[7].split('=')[1])
    
            CGI_Dict[CGIname][5] += 1
            CGI_Dict[CGIname][6] += phastConsBp
            CGI_Dict[CGIname][8] += LOD
            if LOD > CGI_Dict[CGIname][7]:
                CGI_Dict[CGIname][7] = LOD
inPhastConsCGIs.close()

#### WRITE TO OUTPUT
# Peak_Dict[peakName]['basicInfo'] = [lengthPeak_Roller, peakLiftsToHuman, lengthPeak_Human, noFeaturesInHuman, numCpGs_Peak, GC_Peak]
# Peak_Dict[peakName]['CGI'] = [CGIname1, CGIname2, ...]
# Peak_Dict[peakName]['Peak'] = [RPM, RPKM]
# Peak_Dict[peakName]['phastCons'] = [phastCons_number, phastConsBp, phastCons_maxLOD, phastCons_sumLOD]
# Peak_Dict[peakName]['other'] = [OldestSegment_Peak, RepeatBp_Peak]

# CGI_Dict[CGIname] = [lengthInRoller, NumCpGs_Roller, liftsToHuman, lengthInHuman, OldestSegment_CGI,
#                      phastConsCGI_number, phastConsCGI_bp, phastConsCGI_maxLOD, phastConsCGI_sumLOD]

# write header
print('\t'.join(['peakName', 'lengthPeak_Roller', 'peakLiftsToHuman', 'lengthPeak_Human', 'noFeaturesInHuman', 'Peak_numCpGs', 'GC_Peak', 'RPM', 'RPKM', 'Peak_phastConsNumber', 'Peak_phastConsBp', 'Peak_phastConsMaxLOD', 'Peak_phastConsSumLOD', 'Peak_OldestSegment', 'Peak_RepeatBp', 'CGI', 'CGI_number', 'CGI_Bp', 'CGI_numCpGs', 'CGI_numberLiftToHuman', 'CGI_BpInHuman', 'CGI_phastConsNumber', 'CGI_phastConsBp', 'CGI_phastConsMaxLOD', 'CGI_phastConsSum', 'CGI_OldestSegment']))

# print to output for each peakName
for peakName in Peak_Dict:
    outputArray = []
    outputArray.append(peakName)
    
    # append basic info - lengthPeak_Roller, peakLiftsToHuman, lengthPeak_Human, noFeaturesInHuman, numCpGs_Peak
    outputArray += Peak_Dict[peakName]['basicInfo']
    
    # append RPM and RPKM
    outputArray += Peak_Dict[peakName]['Peak']
    
    # append phastCons info
    outputArray += Peak_Dict[peakName]['phastCons']
    
    # append other info
    outputArray += Peak_Dict[peakName]['other']
    
    # append CGI info - requires integrating across potentially multiple CGIs
    numberOfCGIs = len(Peak_Dict[peakName]['CGI'])
    if numberOfCGIs == 0:
        outputArray.append('0')
        outputArray += ['NA'] * 10
    else:
        outputArray.append('1')
        outputArray.append(numberOfCGIs)
        
        CGI_Bp = 0
        CGI_numCpGs = 0
        CGI_numberLiftToHuman = 0
        CGI_BpInHuman = 0
        CGI_phastConsNumber = 0
        CGI_phastConsBp = 0
        CGI_phastConsMaxLOD = 0
        CGI_phastConsSum = 0
        CGI_OldestSegment = 'None'
        
        for CGIname in Peak_Dict[peakName]['CGI']:
            CGI_Bp += CGI_Dict[CGIname][0]
            CGI_numCpGs += CGI_Dict[CGIname][1]
            CGI_numberLiftToHuman += CGI_Dict[CGIname][2]
            CGI_BpInHuman +=  CGI_Dict[CGIname][3]
            CGI_phastConsNumber +=  CGI_Dict[CGIname][5]
            CGI_phastConsBp +=  CGI_Dict[CGIname][6]
            
            if CGI_Dict[CGIname][7] > CGI_phastConsMaxLOD:
                CGI_phastConsMaxLOD = CGI_Dict[CGIname][7]
            
            CGI_phastConsSum +=  CGI_Dict[CGIname][8]
            
            if AgeSegmentation_Index[CGI_Dict[CGIname][4]] > AgeSegmentation_Index[CGI_OldestSegment]:
                CGI_OldestSegment = CGI_Dict[CGIname][4]
    
        outputArray += [CGI_Bp, CGI_numCpGs, CGI_numberLiftToHuman, CGI_BpInHuman]
        outputArray += [CGI_phastConsNumber, CGI_phastConsBp, CGI_phastConsMaxLOD, CGI_phastConsSum]
        outputArray += [CGI_OldestSegment]
        
    outputArrayAsStrings = []
    
    for i in outputArray:
        outputArrayAsStrings.append(str(i))
        
    print('\t'.join(outputArrayAsStrings))


