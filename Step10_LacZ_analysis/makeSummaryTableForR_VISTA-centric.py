# 4/21/23
# Purpose: take intersection files (VISTA with peaks, VISTA with CGIs)
# and output a summary table for use in R
# where each row is a VISTA element
# and columns are: VISTAname, peakNames, peakBinary, LacZ_result, peakLength, peakOverlap, LacZLength, CGI, CGILength, CGI_overlap
# where LacZ_result = sameTissue/otherTissue/negative
# and LacZLength is for the element highest up on the list of sameTissue/otherTissue/negative and if there is a tie, then the biggest one

import sys

inTissue = str(sys.argv[1])
inMark = str(sys.argv[2])

inPeakOverlap = open('intersect_VISTA-centric/e11.5_'+inTissue+'_'+inMark+'_intersectPeaks.txt','rt')
inCGIoverlap = open('intersect_VISTA-centric/e11.5_'+inTissue+'_'+inMark+'_intersectCGIs.txt','rt')


# read in LacZ intersection file and store info in LacZ_Dict

LacZ_Dict = {}
# LacZ_Dict[VISTAname] = [peakNames, peakBinary, LacZ_result, peakLength, peakOverlap, LacZLength]

ResultPriority_Dict = {'sameTissue':3,'otherTissue':2,'negative':1}

for line in inPeakOverlap:
    splitLine = line.strip().split('\t')
    VISTAname = splitLine[3].split(';')[0]
    peakName = splitLine[7]

    # store LacZ_result
    if splitLine[3].split(';')[1] == 'negative':
        LacZ_result = 'negative'
    elif splitLine[3].split(';')[1] == 'positive':
        tissues = splitLine[3].split(';')[2].split(',')
        if inTissue in tissues:
            LacZ_result = 'sameTissue'
        else:
            LacZ_result = 'otherTissue'

    LacZLength = int(splitLine[2]) - int(splitLine[1])
    
    peakBinary = 0
    peakLength = 0
    peakOverlap = 0
    if peakName != '.':
        peakLength = int(splitLine[6]) - int(splitLine[5])
        peakOverlap = int(splitLine[14])
        peakBinary = 1
    
    # store info
    if VISTAname not in LacZ_Dict:
        LacZ_Dict[VISTAname] = [peakName, peakBinary, LacZ_result, peakLength, peakOverlap, LacZLength]
    
    # if info already exists in dictionary, modify info
    elif VISTAname in LacZ_Dict:
        # add new peakName after comma in field 0
        oldPeakName = LacZ_Dict[VISTAname][0]
        LacZ_Dict[VISTAname][0] = oldPeakName + ',' + peakName
        # add to peakLength and peakOverlap
        LacZ_Dict[VISTAname][3] += peakLength
        LacZ_Dict[VISTAname][4] += peakOverlap
inPeakOverlap.close()


# read in CGI intersection file and store info in CGI_Dict (only if peakName in LacZ_Dict)

CGI_Dict = {}
# CGI_Dict[VISTAname] = [CGI, CGILength, CGI_overlap]

for line in inCGIoverlap:
    splitLine = line.strip().split('\t')
    VISTAname = splitLine[3].split(';')[0]
    
    # store info under VISTA element
    CGILength = int(splitLine[6]) - int(splitLine[5])
    CGI_overlap = int(splitLine[7])
    
    # create new entry
    if VISTAname not in CGI_Dict:
        CGI_Dict[VISTAname] = [1, CGILength, CGI_overlap]
    else:
        # or add to existing entry if multiple CGIs
        CGI_Dict[VISTAname][1] += CGILength
        CGI_Dict[VISTAname][2] += CGI_overlap
inCGIoverlap.close()


# go through LacZ_Dict and print to output, incorporating info from CGI_Dict
# VISTAname, peakNames, peakBinary, LacZ_result, peakLength, peakOverlap, LacZLength, CGI, CGILength, CGI_overlap
# LacZ_Dict[VISTAname] = [peakNames, peakBinary, LacZ_result, peakLength, peakOverlap, LacZLength]
# CGI_Dict[VISTAname] = [CGI, CGILength, CGI_overlap]

print('VISTAname'+'\t'+'peakNames'+'\t'+'peakBinary'+'\t'+'LacZ_result'+'\t'+'peakLength'+'\t'+'peakOverlap'+'\t'+'LacZLength'+'\t'+'CGI'+'\t'+'CGILength'+'\t'+'CGI_overlap')

for VISTAname in LacZ_Dict:

    if VISTAname not in CGI_Dict:
        CGI_Dict[VISTAname] = [0,0,0]

    # turn all fields into strings
    outputString = []
    for i in LacZ_Dict[VISTAname]:
        outputString.append(str(i))
    for j in CGI_Dict[VISTAname]:
        outputString.append(str(j))

    print(VISTAname+'\t'+'\t'.join(outputString))
    

