# 3/20/24
# Purpose: take intersection files (cCREs with peaks/CGIs/interactions)
# and output a summary table for use in R
# where each row is a cCRE
# and columns are: cCRE, peakNames, peakBinary, peakOverlap, cCRE_length, CGI, CGILength, CGI_overlap, nameInteractions, numberInteractions
# where LacZ_result = sameTissue/otherTissue/negative
# and LacZLength is for the element highest up on the list of sameTissue/otherTissue/negative and if there is a tie, then the biggest one

import sys

inMark = str(sys.argv[1])

inPeakOverlap = open('intersect_peaks/e14.5_liver_'+inMark+'_intersectPeaks.txt','rt')
inCGIoverlap = open('intersect_CGIs/encodeCcreCombined_mm10_noFeatures_intersectCGIs.txt','rt')
inInteractionOverlap = open('intersect_interactions/encodeCcreCombined_mm10_noFeatures_intersectInteractions.txt')

Peak_Dict = {}
# Peak_Dict[cCRE] = [peakNames, peakBinary, peakOverlap, cCRE_length]

for line in inPeakOverlap:
    splitLine = line.strip().split('\t')
    cCRE = splitLine[3]
    peakName = splitLine[7]
    cCRE_length = int(splitLine[2]) - int(splitLine[1])

    peakBinary = 0
    peakLength = 0
    peakOverlap = 0
    if peakName != '.':
        peakLength = int(splitLine[6]) - int(splitLine[5])
        peakOverlap = int(splitLine[8])
        peakBinary = 1
        
    # store info
    if cCRE not in Peak_Dict:
        Peak_Dict[cCRE] = [peakName, peakBinary, peakOverlap, cCRE_length]
    
    # if info already exists in dictionary, modify info
    elif cCRE in Peak_Dict:
        
        # add new peakName after comma in field 0
        oldPeakName = Peak_Dict[cCRE][0]
        Peak_Dict[cCRE][0] = oldPeakName + ',' + peakName
        
        # add to peakOverlap
        Peak_Dict[cCRE][2] += peakOverlap
inPeakOverlap.close()


# read in CGI intersection file and store info in CGI_Dict

CGI_Dict = {}
# CGI_Dict[cCRE] = [CGI, CGILength, CGI_overlap]

for line in inCGIoverlap:
    splitLine = line.strip().split('\t')
    cCRE = splitLine[3]
    
    # get cCRE info
    CGILength = int(splitLine[6]) - int(splitLine[5])
    CGI_overlap = int(splitLine[7])
    
    # create new entry
    if CGILength > 0:
        if cCRE not in CGI_Dict:
            CGI_Dict[cCRE] = [1, CGILength, CGI_overlap]
        
        else:
            # or add to existing entry if multiple CGIs
            CGI_Dict[cCRE][1] += CGILength
            CGI_Dict[cCRE][2] += CGI_overlap
inCGIoverlap.close()


# read in interactions and store info in Interaction_Dict

Interaction_Dict = {}
# Interaction_Dict[cCRE] = [nameInteractions, numberInteractions]

for line in inInteractionOverlap:
    splitLine = line.strip().split('\t')
    cCRE = splitLine[3]
    nameInteraction = splitLine[7]
    
    # store info
    if nameInteraction != '.':
    
        # new entry
        if cCRE not in Interaction_Dict:
            Interaction_Dict[cCRE] = [nameInteraction, 1]
            
        # or modify existing entry
        else:
            # add new peakName after comma in field 0
            oldInteractionName = Interaction_Dict[cCRE][0]
            Interaction_Dict[cCRE][0] = oldInteractionName + ',' + nameInteraction
        
            # add to number
            Interaction_Dict[cCRE][1] += 1


# go through LacZ_Dict and print to output, incorporating info from CGI_Dict
# cCRE, peakNames, peakBinary, peakOverlap, cCRE_length, CGI, CGILength, CGI_overlap, nameInteractions, numberInteractions
# Peak_Dict[cCRE] = [peakNames, peakBinary, peakOverlap, cCRE_length]
# CGI_Dict[cCRE] = [CGI, CGILength, CGI_overlap]
# Interaction_Dict[cCRE] = [nameInteractions, numberInteractions]

print('cCRE'+'\t'+'peakNames'+'\t'+'peakBinary'+'\t'+'peakOverlap'+'\t'+'cCRE_length'+'\t'+'CGI'+'\t'+'CGILength'+'\t'+'CGI_overlap'+'\t'+'nameInteractions'+'\t'+'numberInteractions')

for cCRE in Peak_Dict:

    if cCRE not in CGI_Dict:
        CGI_Dict[cCRE] = [0,0,0]
        
    if cCRE not in Interaction_Dict:
        Interaction_Dict[cCRE] = [0,0]

    # turn all fields into strings
    outputString = []
    for i in Peak_Dict[cCRE]:
        outputString.append(str(i))
    for j in CGI_Dict[cCRE]:
        outputString.append(str(j))
    for k in Interaction_Dict[cCRE]:
        outputString.append(str(k))

    print(cCRE+'\t'+'\t'.join(outputString))
    

