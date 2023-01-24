# 7/31/22
# Purpose: make table for a given two-species comparison where each row is:
# CGIname CGI_A CGI_B lengthCGI_A lengthCGI_B lengthCGI_human lengthOriginalCGI_A lengthOriginalCGI_B RepeatBasesCGI_A RepeatBasesCGI_B PhastConsBases_CGI MaxPhastConsScore_CGI SumOfPhastConsScores_CGI OldestSegment_CGI GC_CGI_A GC_CGI_B NumCpGs_CGI_A NumCpGs_CGI_B NumCpGs_ReconciledCGI_A NumCpGs_ReconciledCGI_B Peak_A Peak_B PeakLifts PeakLength_Human PhastConsBases_Peak MaxPhastConsScore_Peak SumOfPhastConsScores_Peak OldestSegment_Peak PeakLength_A PeakLength_B GC_Peak_A GC_Peak_B NumCpGs_Peak_A NumCpGs_Peak_B RPM_A RPM_B RPKM_A RPKM_B phastBiasCGI_A phastBiasCGITotalLength_A phastBiasCGIOverlap_A phastBiasCGI_B phastBiasCGITotalLength_B phastBiasCGIOverlap_B phastBiasPeak_A phastBiasPeakTotalLength_A phastBiasPeakOverlap_A phastBiasPeak_B phastBiasPeakTotalLength_B phastBiasPeakOverlap_B
# edited 9/15/22 to include repeat filter and phastBias overlaps
# updated 1/7/22 to include CpG number and length in each original species, rather than from reconciled interval

import sys

inCGI_A = open(sys.argv[1],'rt')               # speciesA_reconciledCGI_intersectSpeciesA_CGI.txt
inCGI_B = open(sys.argv[2],'rt')               # speciesB_reconciledCGI_intersectSpeciesB_CGI.txt
inRepeats_A = open(sys.argv[3],'rt')           # speciesA_reconciledCGI_intersectSpeciesA_rmsk.txt
inRepeats_B = open(sys.argv[4],'rt')           # speciesA_reconciledCGI_intersectSpeciesB_rmsk.txt
inPhastConsCGI = open(sys.argv[5],'rt')        # human_reconciledCGI_intersectPhastCons.txt
inAgeSegmentsCGI = open(sys.argv[6],'rt')      # human_reconciledCGI_intersectAgeSegments.txt
inPeaks_A = open(sys.argv[7],'rt')             # speciesA_reconciledCGI_intersectH3K4me3.txt # this is before getting reconciledPeaks because some may not pass that step
inPeaks_B = open(sys.argv[8],'rt')             # speciesB_reconciledCGI_intersectH3K4me3.txt
inPhastConsPeaks = open(sys.argv[9],'rt')      # human_reconciledPeaks_intersectPhastCons.txt
inAgeSegmentsPeaks = open(sys.argv[10],'rt')   # human_reconciledPeaks_intersectAgeSegments.txt
inReconciledPeaks_A = open(sys.argv[11],'rt')  # speciesA_reconciledCGI_intersect_reconciledPeaks.txt
inReconciledPeaks_B = open(sys.argv[12],'rt')  # speciesB_reconciledCGI_intersect_reconciledPeaks.txt
inSignal_A = open(sys.argv[13],'rt')           # speciesA_H3K4me3.signal
inSignal_B = open(sys.argv[14],'rt')           # speciesB_H3K4me3.signal
inFaCountCGI_A = open(sys.argv[15],'rt')       # speciesA_reconciledCGI.faCount
inFaCountCGI_B = open(sys.argv[16],'rt')       # speciesB_reconciledCGI.faCount
inFaCountPeaks_A = open(sys.argv[17],'rt')     # speciesA_reconciledPeaks.faCount
inFaCountPeaks_B =  open(sys.argv[18],'rt')    # speciesB_reconciledPeaks.faCount

whetherToFilter = str(sys.argv[19]) # repeatFilter vs noRepeatFilter
# filter criteria:
lengthFilter = 1.25 # length in A / length in B < 1.25 (and reciprocal)
repeatFilter = 0.2 # no more than 20% of length is repeat
nonRepeatLengthFilter = 200 # total non-repeat length > 200bp in both species

if len(sys.argv) > 20:
    inPhastBias_CGI_A = open(sys.argv[20],'rt')
    inPhastBias_CGI_B = open(sys.argv[21],'rt')
    inPhastBias_Peaks_A = open(sys.argv[22],'rt')
    inPhastBias_Peaks_B = open(sys.argv[23],'rt')

#### Store reconciled CGI name, length in each species and NumCpGs in each species

ReconciledCGI_Dict = {}
# ReconciledCGI_Dict[name][A/B] = [length, numCpGs from sum of CGIs, # bases covered by repeats, whether peak (0 vs 1),
#                                  whether reconciledPeak (0 vs name), numCpGs from faCount, GC %, length in original genome]
# ReconciledCGI_Dict[name][human] = [length in human, bp covered by phastCons, max phastCons score, sum of phastCons scores, oldest age in segmentation maps]

ReconciledPeaks_Dict = {}
# ReconciledPeaks_Dict[name][A/B] = [length, RPM, RPKM, numCpGs from faCount, GC %]
# ReconciledPeaks_Dict[name][human] = [length in human, bp covered by phastCons, max phastCons score, sum of phastCons scores, oldest age in segmentation maps]

# store info on overlap of reconciled CGIs with phastCons bp

for line in inPhastConsCGI:
    splitLine = line.strip().split()
    name = splitLine[3]
    lengthHuman = int(splitLine[2]) - int(splitLine[1])
    overlap = int(splitLine[8])
    if name not in ReconciledCGI_Dict:
        ReconciledCGI_Dict[name] = {}
        ReconciledCGI_Dict[name]['human'] = [lengthHuman,0,0,0,'None']
        ReconciledCGI_Dict[name]['A'] = [0,0,0,0,0,0,0,0]
        ReconciledCGI_Dict[name]['B'] = [0,0,0,0,0,0,0,0]
    ReconciledCGI_Dict[name]['human'][1] += overlap
    if splitLine[4] != '.':
        lodScore = int(splitLine[7].split('=')[1])
        if lodScore > ReconciledCGI_Dict[name]['human'][2]:
            ReconciledCGI_Dict[name]['human'][2] = lodScore
        ReconciledCGI_Dict[name]['human'][3] += lodScore
inPhastConsCGI.close()

# store info on overlap of reconciled CGIs with human age segmentation map (Emera et al 2016)
OlderThanAmniota = ['Vertebrata','Gnathostomata','Tetrapoda']
AgeOrder = {'OlderThanAmniota':1,'Amniota':2,'Mammalia':3,'Theria':4,'Eutheria':5,'Primate':6,'Ape':7,'Human':8,'None':9}

for line in inAgeSegmentsCGI:
    splitLine = line.strip().split()
    name = splitLine[3]
    overlap = int(splitLine[13])
    age = splitLine[7]
    if age in OlderThanAmniota:
        age = 'OlderThanAmniota'
    if overlap >= 100 and 'NONE' not in age and age != '.':
        if AgeOrder[age] < AgeOrder[ReconciledCGI_Dict[name]['human'][4]]:
            ReconciledCGI_Dict[name]['human'][4] = age
inAgeSegmentsCGI.close()

# move on to other set of files with info in species A and species B

arrayAB = ['A','B']
intersectCGIFiles = [inCGI_A, inCGI_B]
repeatFiles = [inRepeats_A, inRepeats_B]
peakFiles = [inPeaks_A, inPeaks_B]
reconciledPeakFiles = [inReconciledPeaks_A, inReconciledPeaks_B]
faCountCGIFiles = [inFaCountCGI_A, inFaCountCGI_B]

if len(sys.argv) > 20:
    phastBiasCGIFiles = [inPhastBias_CGI_A,inPhastBias_CGI_B]
    phastBiasPeakFiles = [inPhastBias_Peaks_A,inPhastBias_Peaks_B]
    PhastBiasCGI_Dict = {}
    # PhastBiasCGI_Dict[reconciledCGIname]['A'/'B'] = [0/1 for phastBias overlapping CGI, total length of phastBias that overlaps CGI, overlapping bp with CGI]
    PhastBiasPeak_Dict = {}
    # PhastBiasCGI_Dict[reconciledPeakName]['A'/'B'] = [0/1 for phastBias overlapping PEAK, total length of phastBias that overlaps PEAK, overlapping bp with PEAK]

for index in [0,1]:

    # store number of CpGs in CGIs that the reconciled CGI overlaps (may extend beyond boundaries of reconciled CGI)
    # and length of region in each species
    for line in intersectCGIFiles[index]:
        splitLine = line.strip().split('\t')
        name = splitLine[3]
        length = int(splitLine[2]) - int(splitLine[1])
        if splitLine[4] != '.':
            numCpGs = int(splitLine[7].split(': ')[1])
            ReconciledCGI_Dict[name][arrayAB[index]][1] += numCpGs
            lengthOriginal = int(splitLine[6]) - int(splitLine[5])
            ReconciledCGI_Dict[name][arrayAB[index]][7] += lengthOriginal
        else:
            numCpGs = 0
        ReconciledCGI_Dict[name][arrayAB[index]][0] = length
    
    # store info on overlap with RepeatMasker track
    for line in repeatFiles[index]:
        splitLine = line.strip().split('\t')
        reconciledCGIname = splitLine[3]
        overlap = int(splitLine[8])
        ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][2] += overlap
    repeatFiles[index].close()

    # store info on peak presence/absence (not reconciled peaks, just peaks in each genome)
    for line in peakFiles[index]:
        splitLine = line.strip().split('\t')
        reconciledCGIname = splitLine[3]
        if splitLine[4] != '.':
            ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][3] = 1
    peakFiles[index].close()

    # store intersection with reconciledPeaks
    for line in reconciledPeakFiles[index]:
        splitLine = line.strip().split('\t')
        if splitLine[4] != '.':
            reconciledCGIname = splitLine[3]
            reconciledPeakName = splitLine[7]
            
            if ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][4] == 0:
                ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][4] = reconciledPeakName
            else:
                ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][4] += ';'
                ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][4] += reconciledPeakName
                    
    # store info on CpG number and GC% in CGIs
    for line in faCountCGIFiles[index]:
        if line[0] != '#' and line.strip().split('\t')[0] != 'total':
            splitLine = line.strip().split('\t')
            reconciledCGIname = splitLine[0].split('::')[0]
            length = int(splitLine[1])
            numCpGs = int(splitLine[7])
            numGplusC = int(splitLine[3]) + int(splitLine[4])
            GCpercent = round(float(numGplusC)/float(length)*100,2)
            ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][5] = numCpGs
            ReconciledCGI_Dict[reconciledCGIname][arrayAB[index]][6] = GCpercent

    if len(sys.argv) > 20:
        for line in phastBiasCGIFiles[index]:
            splitLine = line.strip().split('\t')
            reconciledCGIname = splitLine[3]
            if reconciledCGIname not in PhastBiasCGI_Dict:
                PhastBiasCGI_Dict[reconciledCGIname] = {}
                PhastBiasCGI_Dict[reconciledCGIname]['A'] = [0,0,0]
                PhastBiasCGI_Dict[reconciledCGIname]['B'] = [0,0,0]
            phastBiasLength = int(splitLine[6]) - int(splitLine[5])
            overlap = int(splitLine[8])

            # add CGI info to PhastBiasCGI_Dict
            if splitLine[4] != '.':
                PhastBiasCGI_Dict[reconciledCGIname][arrayAB[index]][0] = 1
                PhastBiasCGI_Dict[reconciledCGIname][arrayAB[index]][1] += phastBiasLength
                PhastBiasCGI_Dict[reconciledCGIname][arrayAB[index]][2] += overlap
        phastBiasCGIFiles[index].close()

        for line in phastBiasPeakFiles[index]:
            splitLine = line.strip().split('\t')
            reconciledPeakName = splitLine[3]
            if reconciledPeakName not in PhastBiasPeak_Dict:
                PhastBiasPeak_Dict[reconciledPeakName] = {}
                PhastBiasPeak_Dict[reconciledPeakName]['A'] = [0,0,0]
                PhastBiasPeak_Dict[reconciledPeakName]['B'] = [0,0,0]
            phastBiasLength = int(splitLine[6]) - int(splitLine[5])
            overlap = int(splitLine[8])

            # add CGI info to PhastBiasCGI_Dict
            if splitLine[4] != '.':
                PhastBiasPeak_Dict[reconciledPeakName][arrayAB[index]][0] = 1
                PhastBiasPeak_Dict[reconciledPeakName][arrayAB[index]][1] += phastBiasLength
                PhastBiasPeak_Dict[reconciledPeakName][arrayAB[index]][2] += overlap
        phastBiasPeakFiles[index].close()


###########################################################################
# TRANSITION TO STORING INFO ON RECONCILED PEAKS
# already stored intersection of phastBias with reconciled peaks, directly above

# store phastCons bp in human reconciledPeaks
for line in inPhastConsPeaks:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    length = int(splitLine[2]) - int(splitLine[1])
    overlap = int(splitLine[8])
    
    if peakName not in ReconciledPeaks_Dict:
        ReconciledPeaks_Dict[peakName] = {}
        ReconciledPeaks_Dict[peakName]['human'] = [length,0,0,0,'None']
        ReconciledPeaks_Dict[peakName]['A'] = ['.','.','.','.','.']
        ReconciledPeaks_Dict[peakName]['B'] = ['.','.','.','.','.']
    ReconciledPeaks_Dict[peakName]['human'][1] += overlap

    if splitLine[4] != '.':
        lodScore = int(splitLine[7].split('=')[1])
        if lodScore > ReconciledPeaks_Dict[peakName]['human'][2]:
            ReconciledPeaks_Dict[peakName]['human'][2] = lodScore
        ReconciledPeaks_Dict[peakName]['human'][3] += lodScore

# store info on overlap of reconciled peaks with human age segmentation map (Emera et al 2016)
OlderThanAmniota = ['Vertebrata','Gnathostomata','Tetrapoda']
AgeOrder = {'OlderThanAmniota':1,'Amniota':2,'Mammalia':3,'Theria':4,'Eutheria':5,'Primate':6,'Ape':7,'Human':8,'None':9}

for line in inAgeSegmentsPeaks:
    splitLine = line.strip().split()
    peakName = splitLine[3]
    overlap = int(splitLine[13])
    age = splitLine[7]
    if age in OlderThanAmniota:
        age = 'OlderThanAmniota'
    if overlap >= 100 and 'NONE' not in age and age != '.':
        if AgeOrder[age] < AgeOrder[ReconciledPeaks_Dict[peakName]['human'][4]]:
            ReconciledPeaks_Dict[peakName]['human'][4] = age
inAgeSegmentsPeaks.close()
    
signalFiles = [inSignal_A, inSignal_B]
faCountPeaksFiles = [inFaCountPeaks_A, inFaCountPeaks_B]

for index in [0,1]:
    
    # store info on RPM, RPKM
    for line in signalFiles[index]:
        splitLine = line.strip().split('\t')
        peakName = splitLine[0]
        #print(peakName)
        RPM = splitLine[1]
        RPKM = splitLine[2]

        ReconciledPeaks_Dict[peakName][arrayAB[index]][1] = RPM
        ReconciledPeaks_Dict[peakName][arrayAB[index]][2] = RPKM

    # store info on CpG number and GC content
    for line in faCountPeaksFiles[index]:
        if line[0] != '#' and line.strip().split('\t')[0] != 'total':
            splitLine = line.strip().split('\t')
            peakName = splitLine[0].split('::')[0]
            length = int(splitLine[1])
            numCpGs = int(splitLine[7])
            numGplusC = int(splitLine[3]) + int(splitLine[4])
            GCpercent = round(float(numGplusC)/float(length)*100,2)
            ReconciledPeaks_Dict[peakName][arrayAB[index]][0] = length
            ReconciledPeaks_Dict[peakName][arrayAB[index]][3] = numCpGs
            ReconciledPeaks_Dict[peakName][arrayAB[index]][4] = GCpercent

#############################################################################
# COLLECT INFO AND WRITE TO OUTPUT
# CGIname CGI_A CGI_B lengthCGI_A lengthCGI_B lengthCGI_human lengthOriginalCGI_A lengthOriginalCGI_B RepeatBasesCGI_A RepeatBasesCGI_B PhastConsBases_CGI MaxPhastConsScore_CGI SumOfPhastConsScores_CGI OldestSegment_CGI GC_CGI_A GC_CGI_B NumCpGs_CGI_A NumCpGs_CGI_B NumCpGs_ReconciledCGI_A NumCpGs_ReconciledCGI_B Peak_A Peak_B PeakLifts PeakLength_Human PhastConsBases_Peak MaxPhastConsScore_Peak SumOfPhastConsScores_Peak OldestSegment_Peak PeakLength_A PeakLength_B GC_Peak_A GC_Peak_B NumCpGs_Peak_A NumCpGs_Peak_B RPM_A RPM_B RPKM_A RPKM_B phastBiasCGI_A phastBiasCGITotalLength_A phastBiasCGIOverlap_A phastBiasCGI_B phastBiasCGITotalLength_B phastBiasCGIOverlap_B phastBiasPeak_A phastBiasPeakTotalLength_A phastBiasPeakOverlap_A phastBiasPeak_B phastBiasPeakTotalLength_B phastBiasPeakOverlap_B

headerArray = ['CGIname','CGI_A','CGI_B','lengthCGI_A','lengthCGI_B','lengthCGI_human','lengthOriginalCGI_A','lengthOriginalCGI_B','RepeatBasesCGI_A','RepeatBasesCGI_B','PhastConsBases_CGI','MaxPhastConsScore_CGI','SumOfPhastConsScores_CGI','OldestSegment_CGI','GC_CGI_A','GC_CGI_B','NumCpGs_CGI_A','NumCpGs_CGI_B','NumCpGs_ReconciledCGI_A','NumCpGs_ReconciledCGI_B','Peak_A','Peak_B','PeakLifts','PeakLength_Human','PhastConsBases_Peak','MaxPhastConsScore_Peak','SumOfPhastConsScores_Peak','OldestSegment_Peak','PeakLength_A','PeakLength_B','GC_Peak_A','GC_Peak_B','NumCpGs_Peak_A','NumCpGs_Peak_B','RPM_A','RPM_B','RPKM_A','RPKM_B','phastBiasCGI_A','phastBiasCGITotalLength_A','phastBiasCGIOverlap_A','phastBiasCGI_B','phastBiasCGITotalLength_B','phastBiasCGIOverlap_B','phastBiasPeak_A','phastBiasPeakTotalLength_A','phastBiasPeakOverlap_A','phastBiasPeak_B','phastBiasPeakTotalLength_B','phastBiasPeakOverlap_B']

if len(sys.argv) > 20:
    print('\t'.join(headerArray))
else:
    print('\t'.join(headerArray[:-12]))

for reconciledCGIname in ReconciledCGI_Dict:

    # FILTER BASED ON LENGTH AND REPEAT CONTENT
    meetsFinalRepeatFilter = 'yes'
    
    if whetherToFilter == 'repeatFilter':
        ratioAoverB = float(ReconciledCGI_Dict[reconciledCGIname]['A'][0]) / float(ReconciledCGI_Dict[reconciledCGIname]['B'][0])
        ratioBoverA = float(ReconciledCGI_Dict[reconciledCGIname]['B'][0]) / float(ReconciledCGI_Dict[reconciledCGIname]['A'][0])
        repeatPercentA = float(ReconciledCGI_Dict[reconciledCGIname]['A'][2]) / float(ReconciledCGI_Dict[reconciledCGIname]['A'][0])
        repeatPercentB = float(ReconciledCGI_Dict[reconciledCGIname]['B'][2]) / float(ReconciledCGI_Dict[reconciledCGIname]['B'][0])
        nonRepeatLengthA = ReconciledCGI_Dict[reconciledCGIname]['A'][0] - ReconciledCGI_Dict[reconciledCGIname]['A'][2]
        nonRepeatLengthB = ReconciledCGI_Dict[reconciledCGIname]['B'][0] - ReconciledCGI_Dict[reconciledCGIname]['B'][2]

        if ratioAoverB > lengthFilter or ratioBoverA > lengthFilter or repeatPercentA > repeatFilter or repeatPercentB > repeatFilter or nonRepeatLengthA < nonRepeatLengthFilter or nonRepeatLengthB < nonRepeatLengthFilter:
            meetsFinalRepeatFilter = 'no'

    if meetsFinalRepeatFilter == 'yes':
            
        outputString = []

        # CGIname, CGI_A, CGI_B
        outputString.append(reconciledCGIname)
        if 'a_' in reconciledCGIname:
            outputString.append('1')
        else:
            outputString.append('0')
        if 'b_' in reconciledCGIname:
            outputString.append('1')
        else:
            outputString.append('0')

        # lengthCGI_A, lengthCGI_B, lengthCGI_human, lengthOriginalCGI_A lengthOriginalCGI_B
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][0])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][0])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['human'][0])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][7])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][7])

        # RepeatBasesCGI_A, RepeatBasesCGI_B
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][2])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][2])

        # PhastConsBases_CGI, MaxPhastConsScore_CGI, SumOfPhastConsScores_CGI, OldestSegment_CGI
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['human'][1])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['human'][2])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['human'][3])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['human'][4])

        # GC_CGI_A, GC_CGI_B, NumCpGs_CGI_A, NumCpGs_CGI_B, NumCpGs_ReconciledCGI_A, NumCpGs_ReconciledCGI_B
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][6])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][6])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][1])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][1])   
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][5])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][5])

        # Peak_A, Peak_B
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['A'][3])
        outputString.append(ReconciledCGI_Dict[reconciledCGIname]['B'][3])

        #print(ReconciledCGI_Dict[reconciledCGIname]['A'][4])

        # include reconciledPeak in analysis only if a reconciled CGI intersects a single reconciledPeak, and if it's the same reconciledPeak in A and B
        if ReconciledCGI_Dict[reconciledCGIname]['A'][4] != 0 and ';' not in ReconciledCGI_Dict[reconciledCGIname]['A'][4] and ReconciledCGI_Dict[reconciledCGIname]['A'][4] == ReconciledCGI_Dict[reconciledCGIname]['B'][4]:
            reconciledPeakName = ReconciledCGI_Dict[reconciledCGIname]['A'][4]

            # PeakLifts, PeakLength_Human, PhastConsBases_Peak, MaxPhastConsScore_Peak, SumOfPhastConsScores_Peak, OldestSegment_Peak
            outputString += ['1']
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['human'][0])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['human'][1])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['human'][2])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['human'][3])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['human'][4])

            # PeakLength_A, PeakLength_B, GC_Peak_A, GC_Peak_B, NumCpGs_Peak_A, NumCpGs_Peak_B
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['A'][0])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['B'][0])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['A'][4])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['B'][4])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['A'][3])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['B'][3])

            # RPM_A, RPM_B, RPKM_A, RPKM_B
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['A'][1])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['B'][1])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['A'][2])
            outputString.append(ReconciledPeaks_Dict[reconciledPeakName]['B'][2])

        else:
            outputString += ['0']
            outputString += ['NA'] * 15

        # add phastBias info
        if len(sys.argv) > 20:
            # phastBiasCGI_A, phastBiasCGITotalLength_A, phastBiasCGIOverlap_A, phastBiasCGI_B, phastBiasCGITotalLength_B, phastBiasCGIOverlap_B
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['A'][0])
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['A'][1])
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['A'][2])
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['B'][0])
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['B'][1])
            outputString.append(PhastBiasCGI_Dict[reconciledCGIname]['B'][2])

            # phastBiasPeak_A, phastBiasPeakTotalLength_A, phastBiasPeakOverlap_A, phastBiasPeak_B, phastBiasPeakTotalLength_B, phastBiasPeakOverlap_B
            if ReconciledCGI_Dict[reconciledCGIname]['A'][4] != 0 and ';' not in ReconciledCGI_Dict[reconciledCGIname]['A'][4] and ReconciledCGI_Dict[reconciledCGIname]['A'][4] == ReconciledCGI_Dict[reconciledCGIname]['B'][4]:
                reconciledPeakName = ReconciledCGI_Dict[reconciledCGIname]['A'][4]

                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['A'][0])
                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['A'][1])
                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['A'][2])
                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['B'][0])
                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['B'][1])
                outputString.append(PhastBiasPeak_Dict[reconciledPeakName]['B'][2])

            else:
                outputString += ['NA'] * 6

        # convert all entries to strings
        outputStringAsStrings = []
        for i in outputString:
            outputStringAsStrings.append(str(i))

        print('\t'.join(outputStringAsStrings))
