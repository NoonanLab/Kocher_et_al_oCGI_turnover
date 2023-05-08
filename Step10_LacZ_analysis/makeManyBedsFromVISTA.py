# 7/17/22
# Purpose: take VISTA file (from Len Pennacchio) and output bed files
# which are specific to each time point x backbone

# Usage:
# python makeManyBedsFromVISTA.py 2022-07-15_Acadia.tsv tissue.dictionary.csv

import sys

inSpreadsheet = open(sys.argv[1],'rt')
inTranslate = open(sys.argv[2],'rt')

# store info from tissue.dictionary.csv in dictionary to translate tissues later
Tissue_Dict = {}
# Tissue_Dict[tissueAbbreviation] = tissueShort

for line in inTranslate:
    splitLine = line.strip().split(',')
    if splitLine[0] != 'tissue':
        Tissue_Dict[splitLine[-1]] = splitLine[-2]
inTranslate.close()

#for i in Tissue_Dict:
#    print(i+'\t'+Tissue_Dict[i])

# store info from full spreadsheet of all elements in a series of nested dictionaries
Element_Dict = {}
# Element_Dict[element][backbone][timePoint] = [pos_neg,assembly,coord,tissues]

uniqueAssembly = []
uniqueTimePoint = []
uniqueBackbone = []

for line in inSpreadsheet:
    splitLine = line.strip().split('\t')
    #print(splitLine)
    trimmedSplitLine = []
    for i in splitLine:
        trimmedSplitLine.append(i.strip('\"'))
    #print(trimmedSplitLine)

    if trimmedSplitLine[0] != 'vista.id':
        element = trimmedSplitLine[0]
        backbone = trimmedSplitLine[1]
        timePoint = trimmedSplitLine[2]
        pos_neg = trimmedSplitLine[3]
        assembly = trimmedSplitLine[7]
        coord = trimmedSplitLine[8]

        # Tissue_Dict[tissueAbbreviation] = tissueShort
        tissuesAbbrev = trimmedSplitLine[4].split(';')
        tissuesShort = []
        if pos_neg == 'positive':
            for i in tissuesAbbrev:
                tissuesShort.append(Tissue_Dict[i])
            tissues = ','.join(tissuesShort)
        else:
            tissues = '.'

        if assembly not in uniqueAssembly:
            uniqueAssembly.append(assembly)
        if timePoint not in uniqueTimePoint:
            uniqueTimePoint.append(timePoint)
        if backbone not in uniqueBackbone:
            uniqueBackbone.append(backbone)

        if element not in Element_Dict:
            Element_Dict[element] = {}
        if backbone not in Element_Dict[element]:
            Element_Dict[element][backbone] = {}
        if timePoint not in Element_Dict[element][backbone]:
            Element_Dict[element][backbone][timePoint] = [pos_neg,assembly,coord,tissues]
    
inSpreadsheet.close()

#print(uniqueAssembly)
#print(uniqueTimePoint)
#print(uniqueBackbone)

# open an output bed file for each species x timePoint x backbone
# and write bed info to each output file
# column 4 has elementName;pos_neg;tissue1,tissue2

for assembly in uniqueAssembly:
    for timePoint in uniqueTimePoint:
        for backbone in uniqueBackbone:
            outBed = open(assembly+'_'+timePoint+'_'+backbone+'.bed','wt')

            for element in Element_Dict:
                if backbone in Element_Dict[element]:
                    if timePoint in Element_Dict[element][backbone]:
                        if Element_Dict[element][backbone][timePoint][1] == assembly:
                    
                            chrom = Element_Dict[element][backbone][timePoint][2].split(':')[0]
                            start = Element_Dict[element][backbone][timePoint][2].split(':')[1].split('-')[0]
                            end = Element_Dict[element][backbone][timePoint][2].split('-')[1]

                            pos_neg = Element_Dict[element][backbone][timePoint][0]
                            tissues = Element_Dict[element][backbone][timePoint][3]

                            outBed.write(chrom+'\t'+start+'\t'+end+'\t'+element+';'+pos_neg+';'+tissues+'\n')
                    
            outBed.close()
