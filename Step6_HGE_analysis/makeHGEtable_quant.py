# 1/29/23
# Purpose: make table for use in R to understand relationship between CGIs and HGEs
# Modified from makeHGEtable.py - this time incorporate more info on each site

# Intended output:
# PeakName_human PeakCoord_human Enhancer_AK Promoter_AK HGE HGP PeakLength_human CGI_human CGILength_human CpGNum_humanCGI LiftsToRhesus IntervalLength_rhesus CGI_rhesus CGILength_rhesus CpGNum_rhesusCGI LiftsToMouse IntervalLength_mouse CGI_mouse CGILength_mouse CpGNum_mouseCGI CpGNum_humanPeak CpGNum_rhesusInterval CpGNum_mouseInterval avgSignal_rep1 totalSignal_rep1 avgSignal_rep2 totalSignal_rep1

# python makeHGEtable_quant.py <tissue> <mark> <timePoint> > HGEtable_tissue_mark_timePoint.txt

import sys

inTissue = str(sys.argv[1])
inMark = str(sys.argv[2])
inTime = str(sys.argv[3])

if inTissue == 'brain':
    if inMark == 'ac':
        inHumanPeaks = open('/home/ak2267/project/EnhancerClasses/hg19/'+inMark+'/merge_'+inTime+'_overlap_named.bed','rt')
        inHumanPeaksCGIs = open('brain/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectCGIs.txt','rt')
    if inMark == 'me2':
        inHumanPeaks = open('/home/ak2267/project/EnhancerClasses/hg19/'+inMark+'/merge_'+inTime+'_me2_overlap_named.bed','rt')
        inHumanPeaksCGIs = open('brain/hg19/'+inMark+'/merge_'+inTime+'_me2_overlap_intersectCGIs.txt','rt')
if inTissue == 'limb':
    inHumanPeaks = open('/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_'+inTime+'_overlap_named.bed','rt')
    inHumanPeaksCGIs = open('limb/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectCGIs.txt','rt')

# these file paths are relative to running the script from within /home/ak2267/project/EnhancerClasses/HGE
inHumanPeaksPromoters = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectProteinCodingPromoters.bed','rt')
inHumanPeaksIntronicIntergenic = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectIntronicIntergenic.bed','rt')

if inTissue == 'brain':
    inHumanPeaksHGEs = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectHGEs.txt', 'rt')
    inHumanPeaksHGPs = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectHGPs.txt', 'rt')
# only one gain file for limb - make it both HGE and HGP and deal it with later in R
if inTissue == 'limb':
    inHumanPeaksHGEs = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectGains.txt', 'rt')
    inHumanPeaksHGPs = open(inTissue+'/hg19/'+inMark+'/merge_'+inTime+'_overlap_intersectGains.txt', 'rt')

inLiftToRhesus = open(inTissue+'/hg19/'+inMark+'/hg19_'+inMark+'_'+inTime+'_sitesThatMapToRheMac2_inRheMac2.bed','rt')
inLiftToMouse = open(inTissue+'/hg19/'+inMark+'/hg19_'+inMark+'_'+inTime+'_sitesThatMapToMm9_inMm9.bed','rt')

inLiftToRhesusCGIs = open(inTissue+'/hg19/'+inMark+'/hg19_'+inMark+'_'+inTime+'_sitesThatMapToRheMac2_intersectRheMac2CGIs.bed','rt')
inLiftToMouseCGIs = open(inTissue+'/hg19/'+inMark+'/hg19_'+inMark+'_'+inTime+'_sitesThatMapToMm9_intersectMm9CGIs.bed','rt')

inQuantRep1 = open('../quant/'+inTissue+'/hg19/'+inMark+'/'+inTime+'_'+inMark+'_quant_rep1.tab','rt')
inQuantRep2 = open('../quant/'+inTissue+'/hg19/'+inMark+'/'+inTime+'_'+inMark+'_quant_rep2.tab','rt')

inHumanFaCount = open('CpG/'+inTissue+'/hg19/'+inMark+'/'+inTime+'.faCount','rt')
inRhesusFaCount = open('CpG/'+inTissue+'/rheMac2/'+inMark+'/'+inTime+'.faCount','rt')
inMouseFaCount = open('CpG/'+inTissue+'/mm9/'+inMark+'/'+inTime+'.faCount','rt')

# Dictionaries to translate between mergedPeak coordinates and line numbers (both ways), plus store time points where active in human
HumanPeak_Dict = {}
# HumanPeak_Dict[coord] = peakNumber
PeakNumber_Dict = {}
# PeakNumber_Dict[peakNumber] = coord

for line in inHumanPeaks:
    splitLine = line.strip().split('\t')
    coord = splitLine[0] + ':' + splitLine[1] + '-' + splitLine[2]
    peakNumber = splitLine[3]
    HumanPeak_Dict[coord] = peakNumber
    PeakNumber_Dict[peakNumber] = coord
inHumanPeaks.close()
    
# Lists to store peak names for peaks that are in promoters and enhancers
PromoterPeaks = []
# PromoterPeaks = [peakName1,peakName2,...]

EnhancerPeaks = []
# EnhancerPeaks = [peakName1,peakName2,...]

for line in inHumanPeaksPromoters:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    PromoterPeaks.append(peakName)
inHumanPeaksPromoters.close()

for line in inHumanPeaksIntronicIntergenic:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    EnhancerPeaks.append(peakName)
inHumanPeaksIntronicIntergenic.close()

# Dictionary to store CGI info for peaks that have a CGI
HumanCGI_Dict = {}
# HumanCGI_Dict[peakName] = [CGI_length,CpGNum]
for line in inHumanPeaksCGIs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        if peakName not in HumanCGI_Dict:
            HumanCGI_Dict[peakName] = [0,0]
        CGI_length = int(splitLine[6]) - int(splitLine[5])
        CpGNum = int(splitLine[7].split(': ')[1])
        HumanCGI_Dict[peakName][0] += CGI_length
        HumanCGI_Dict[peakName][1] += CpGNum
inHumanPeaksCGIs.close()
    
# Dictionaries to store info for peaks that are HGEs or HGPs
HGE_List = []
# HGE_List = [peakName1, peakName2, ...]
HGP_List = []
# HGP_List = [peakName1, peakName2, ...]

for line in inHumanPeaksHGEs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        if peakName not in HGE_List:
            HGE_List.append(peakName)
inHumanPeaksHGEs.close()
    
for line in inHumanPeaksHGPs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        if peakName not in HGP_List:
            HGP_List.append(peakName)
inHumanPeaksHGPs.close()

# Dictionaries to store peaks that lift to Rhesus and Mouse, and the lengths of those sequences in those species
RhesusLength_Dict = {}
# RhesusLength_Dict[peakName] = length of sequence in rhesus
MouseLength_Dict = {}
# MouseLength_Dict[peakName] = length of sequence in mouse

for line in inLiftToRhesus:
    splitLine = line.strip().split('\t')
    length = int(splitLine[2]) - int(splitLine[1])
    peakName = splitLine[3]
    RhesusLength_Dict[peakName] = length
inLiftToRhesus.close()

for line in inLiftToMouse:
    splitLine = line.strip().split('\t')
    length = int(splitLine[2]) - int(splitLine[1])
    peakName = splitLine[3]
    MouseLength_Dict[peakName] = length
inLiftToMouse.close()

# Dictionaries to store CGI and peak status in Rhesus and Mouse
#inLiftToRhesusCGIs = open('hg19/'+inMark+'/merge_'+inTime+'_overlap_LOrheMac2_intersectRheMac2CGIs.bed','rt')
#inLiftToMouseCGIs = open('hg19/'+inMark+'/merge_'+inTime+'_overlap_LOmm9_intersectMm9CGIs.bed','rt')

RhesusCGI_Dict = {}
# RhesusCGI_Dict[peakName] = [CGI_length,CpGNum]

MouseCGI_Dict = {}
# MouseCGI_Dict[peakName] = [CGI_length,CpGNum]

for line in inLiftToRhesusCGIs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        CGI_length = int(splitLine[6]) - int(splitLine[5])
        CpGNum = int(splitLine[7].split(': ')[1])
        if peakName not in RhesusCGI_Dict:
            RhesusCGI_Dict[peakName] = [0,0]
        RhesusCGI_Dict[peakName][0] += CGI_length
        RhesusCGI_Dict[peakName][1] += CpGNum
inLiftToRhesusCGIs.close()

for line in inLiftToMouseCGIs:
    splitLine = line.strip().split('\t')
    peakName = splitLine[3]
    if splitLine[4] != '.':
        CGI_length = int(splitLine[6]) - int(splitLine[5])
        CpGNum = int(splitLine[7].split(': ')[1])
        if peakName not in MouseCGI_Dict:
            MouseCGI_Dict[peakName] = [0,0]
        MouseCGI_Dict[peakName][0] += CGI_length
        MouseCGI_Dict[peakName][1] += CpGNum
inLiftToMouseCGIs.close()

#inQuantRep1 = open('quant/'+inMark+'/'+inTime+'_'+inMark+'_quant_rep1.tab','rt')
#inQuantRep2 = open('quant/'+inMark+'/'+inTime+'_'+inMark+'_quant_rep2.tab','rt')

# store signal collected by bigWigAverageOverBed
QuantRep1_Dict = {}
# QuantRep1_Dict[peakName] = [sum of signal, average signal]
QuantRep2_Dict = {}
# QuantRep2_Dict[peakName] = [sum of signal, average signal]

for line in inQuantRep1:
    splitLine = line.strip().split('\t')
    peakName = splitLine[0]
    total = splitLine[3]
    mean = splitLine[4]
    QuantRep1_Dict[peakName] = [total,mean]
inQuantRep1.close()

for line in inQuantRep2:
    splitLine = line.strip().split('\t')
    peakName = splitLine[0]
    total = splitLine[3]
    mean = splitLine[4]
    QuantRep2_Dict[peakName] = [total,mean]
inQuantRep2.close()

#inHumanFaCount = open('CpG/hg19/'+inMark+'/'+inTime+'.faCount','rt')
#inRhesusFaCount = open('CpG/rheMac2/'+inMark+'/'+inTime+'.faCount','rt')
#inMouseFaCount = open('CpG/mm9/'+inMark+'/'+inTime+'.faCount','rt')

# store CpG numbers across entire peak interval in each species
PeakCpGs_Dict = {}
PeakCpGs_Dict['Human'] = {}
# PeakCpGs_Dict['Human'][peakName] = # CpGs
PeakCpGs_Dict['Rhesus'] = {}
# PeakCpGs_Dict['Rhesus'][peakName] = # CpGs
PeakCpGs_Dict['Mouse'] = {}
# PeakCpGs_Dict['Mouse'][peakName] = # CpGs

for line in inHumanFaCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        peakName = splitLine[0]
        CpGs = splitLine[7]
        PeakCpGs_Dict['Human'][peakName] = CpGs
inHumanFaCount.close()

for line in inRhesusFaCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        peakName = splitLine[0]
        CpGs = splitLine[7]
        PeakCpGs_Dict['Rhesus'][peakName] = CpGs
inRhesusFaCount.close()

for line in inMouseFaCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        peakName = splitLine[0]
        CpGs = splitLine[7]
        PeakCpGs_Dict['Mouse'][peakName] = CpGs
inMouseFaCount.close()
        
# PeakName_human PeakCoord_human Enhancer_AK Promoter_AK HGE HGP PeakLength_human CGI_human CGILength_human CpGNum_humanCGI LiftsToRhesus IntervalLength_rhesus CGI_rhesus CGILength_rhesus CpGNum_rhesusCGI LiftsToMouse IntervalLength_mouse CGI_mouse CGILength_mouse CpGNum_mouseCGI CpGNum_humanPeak CpGNum_rhesusInterval CpGNum_mouseInterval avgSignal_rep1 totalSignal_rep1 avgSignal_rep2 totalSignal_rep1

# HumanPeak_Dict[coord] = peakNumber
# PeakNumber_Dict[peakNumber] = coord

# PromoterPeaks = [peakName1,peakName2,...]
# EnhancerPeaks = [peakName1,peakName2,...]

# HumanCGI_Dict[peakName] = [CGI_length,CpGNum]
# HGE_List = [peakName1, peakName2, ...]
# HGP_List = [peakName1, peakName2, ...]

# RhesusLength_Dict[peakName] = length of sequence in rhesus
# MouseLength_Dict[peakName] = length of sequence in mouse
# RhesusCGI_Dict[peakName] = [CGI_length,CpGNum]
# MouseCGI_Dict[peakName] = [CGI_length,CpGNum]

# QuantRep1_Dict[peakName] = [sum of signal, average signal]
# QuantRep2_Dict[peakName] = [sum of signal, average signal]
# PeakCpGs_Dict['Human'][peakName] = # CpGs
# PeakCpGs_Dict['Rhesus'][peakName] = # CpGs
# PeakCpGs_Dict['Mouse'][peakName] = # CpGs

print('PeakName_human'+'\t'+'PeakCoord_human'+'\t'+'Enhancer_AK'+'\t'+'Promoter_AK'+'\t'+'HGE'+'\t'+'HGP'+'\t'+'PeakLength_human'+'\t'+'CGI_human'+'\t'+'CGILength_human'+'\t'+'CpGNum_humanCGI'+'\t'+'LiftsToRhesus'+'\t'+'IntervalLength_rhesus'+'\t'+'CGI_rhesus'+'\t'+'CGILength_rhesus'+'\t'+'CpGNum_rhesusCGI'+'\t'+'LiftsToMouse'+'\t'+'IntervalLength_mouse'+'\t'+'CGI_mouse'+'\t'+'CGILength_mouse'+'\t'+'CpGNum_mouseCGI'+'\t'+'CpGNum_humanPeak'+'\t'+'CpGNum_rhesusInterval'+'\t'+'CpGNum_mouseInterval'+'\t'+'avgSignal_rep1'+'\t'+'totalSignal_rep1'+'\t'+'avgSignal_rep2'+'\t'+'totalSignal_rep2')

for peakName in PeakNumber_Dict:
    outputArray = []

    # name and coord
    outputArray.append(peakName)
    coord = PeakNumber_Dict[peakName]
    outputArray.append(coord)

    # whether enhancer or promoter
    if peakName in EnhancerPeaks:
        outputArray.append('1')
    else:
        outputArray.append('0')
    if peakName in PromoterPeaks:
        outputArray.append('1')
    else:
        outputArray.append('0')

    # gain status
    if peakName in HGE_List:
        outputArray.append('1')
    else:
        outputArray.append('0')

    if peakName in HGP_List:
        outputArray.append('1')
    else:
        outputArray.append('0')

    # peak length in human
    length = int(coord.split(':')[1].split('-')[1]) - int(coord.split(':')[1].split('-')[0])
    outputArray.append(str(length))
    
    # CGI status in human, plus CGI length and CpGNum
    if peakName in HumanCGI_Dict:
        outputArray.append('1')
        outputArray.append(str(HumanCGI_Dict[peakName][0]))
        outputArray.append(str(HumanCGI_Dict[peakName][1]))
    else:
        outputArray.append('0')
        outputArray.append('NA')
        outputArray.append('NA')

    # whether it lifts to rhesus and mouse, then length of interval, and length of overlapping peaks
    # and CGI info for rhesus and mouse since it's easy to add here
    # RhesusLength_Dict[peakName] = length of sequence in rhesus
    # MouseLength_Dict[peakName] = length of sequence in mouse
    if peakName in RhesusLength_Dict:
        outputArray.append('1')
        outputArray.append(str(RhesusLength_Dict[peakName]))
        
        # CGI info
        if peakName in RhesusCGI_Dict:
            outputArray.append('1')
            outputArray.append(str(RhesusCGI_Dict[peakName][0]))
            outputArray.append(str(RhesusCGI_Dict[peakName][1]))
        else:
            outputArray.append('0')
            outputArray.append('NA')
            outputArray.append('NA')
    else:
        outputArray.append('0')
        outputArray.append('NA')
        outputArray.append('NA')
        outputArray.append('NA')
        outputArray.append('NA')
        
    if peakName in MouseLength_Dict:
        outputArray.append('1')
        outputArray.append(str(MouseLength_Dict[peakName]))

        # CGI info
        if peakName in MouseCGI_Dict:
            outputArray.append('1')
            outputArray.append(str(MouseCGI_Dict[peakName][0]))
            outputArray.append(str(MouseCGI_Dict[peakName][1]))
        else:
            outputArray.append('0')
            outputArray.append('NA')
            outputArray.append('NA')
    else:
        outputArray.append('0')
        outputArray.append('NA')
        outputArray.append('NA')
        outputArray.append('NA')
        outputArray.append('NA')
    
    # CpG numbers in total peak interval
    outputArray.append(PeakCpGs_Dict['Human'][peakName])
    if peakName in PeakCpGs_Dict['Rhesus']:
        outputArray.append(PeakCpGs_Dict['Rhesus'][peakName])
    else:
        outputArray.append('NA')
    if peakName in PeakCpGs_Dict['Mouse']:
        outputArray.append(PeakCpGs_Dict['Mouse'][peakName])
    else:
        outputArray.append('NA')
        
    # peak signal
    outputArray.append(QuantRep1_Dict[peakName][1])
    outputArray.append(QuantRep1_Dict[peakName][0])
    outputArray.append(QuantRep2_Dict[peakName][1])
    outputArray.append(QuantRep2_Dict[peakName][0])

    print('\t'.join(outputArray))

    
# PeakName_human PeakCoord_human Enhancer_AK Promoter_AK HGE HGP PeakLength_human CGI_human CGILength_human CpGNum_humanCGI LiftsToRhesus IntervalLength_rhesus CGI_rhesus CGILength_rhesus CpGNum_rhesusCGI LiftsToMouse IntervalLength_mouse CGI_mouse CGILength_mouse CpGNum_mouseCGI CpGNum_humanPeak CpGNum_rhesusInterval CpGNum_mouseInterval avgSignal_rep1 totalSignal_rep1 avgSignal_rep2 totalSignal_rep1

# QuantRep1_Dict[peakName] = [sum of signal, average signal]
# QuantRep2_Dict[peakName] = [sum of signal, average signal]
# PeakCpGs_Dict['Human'][peakName] = # CpGs
# PeakCpGs_Dict['Rhesus'][peakName] = # CpGs
# PeakCpGs_Dict['Mouse'][peakName] = # CpGs
