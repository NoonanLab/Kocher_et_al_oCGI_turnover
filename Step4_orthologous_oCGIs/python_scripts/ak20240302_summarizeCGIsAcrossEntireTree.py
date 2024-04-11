# 3/2/2024
# Purpose: take bed files for CGIs in each species and summarize their status in all species in the tree
# Status can be 0 = sequence not present, 1 = sequence present but no CGI, 2 = sequence present and CGI present

# python summarizeCGIsAcrossEntireTree.py rheMac10,calJac4,mm39,rn7,susScr11,canFam6,felCat9,equCab3

import sys

inSpeciesList = str(sys.argv[1]).split(',')

# store species position in list (0-based)
SpeciesIndex_Dict = {}
# SpeciesIndex_Dict[species] = #

for i in range(0,len(inSpeciesList)):
    SpeciesIndex_Dict[inSpeciesList[i]] = i

#for j in SpeciesIndex_Dict:
#    print(str(j)+' '+str(SpeciesIndex_Dict[j]))

# file name ends
fileEnd_liftsToOtherSpecies = '_CGIsThatMapBackToHuman_inOwnCoord.bed'
fileEnd_liftsToOtherSpecies_andNoFeature = '_mergedCGIs_noFeatures_in' # species.bed

# open all CGIs in human
inAllCGIs_human = open('mergedCGIs_hg38.bed', 'rt')

# store which chrom each CGI is on in human for making dictionaries more efficient
chromOfCGI_Dict = {}
# chromOfCGI_Dict[CGInameString] = chrom_in_hg38

# store list of all chrom possibilities
chromList = []

for line in inAllCGIs_human:
    splitLine = line.strip().split('\t')
    chrom_in_hg38 = splitLine[0]
    CGInameString = splitLine[3]
    
    # store chromosome of each CGInameString
    chromOfCGI_Dict[CGInameString] = chrom_in_hg38
    
    # add to running list of all chrom possibilities
    if chrom_in_hg38 not in chromList:
        chromList.append(chrom_in_hg38)
        
print(chromList)

# open phastCons file
inPhastCons = open('mergedCGIs_hg38_overlapPhastCons.bed', 'rt')

## create dictionary all CGInameStrings that have a phastCons element

phastCons_Dict = {}
# phastCons_Dict[chrom]= [CGI1, CGI2, CGI3, ...]

# store chrom possibilities
for chrom in chromList:
    phastCons_Dict[chrom] = []

for line in inPhastCons:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    CGInameString = splitLine[3]
    
    # add CGInameString to phastCons dict to indicate it overlaps a phastCons element
    phastCons_Dict[chrom].append(CGInameString)
    
inPhastCons.close()

# dictionary containing CGIs that lift to each species
CGIs_lift_Dict = {}
# CGIs_lift_Dict[species][chrom_in_hg38] = [CGI1, CGI2, CGI3, ...]

# dictionary containing CGIs that lift to each species AND don't overlap a feature
CGIs_lift_noFeature_Dict = {}
# CGIs_lift_noFeature_Dict[species][chrom_in_hg38] = [CGI1, CGI2, CGI3, ...]

# load both of the above dictionaries

for species in inSpeciesList:

    print(species)
    
    # initialize dictionaries and create nested dictionaries with chromosomes
    CGIs_lift_Dict[species] = {}
    CGIs_lift_noFeature_Dict[species] = {}
    
    for chrom in chromList:
        CGIs_lift_Dict[species][chrom] = []
        CGIs_lift_noFeature_Dict[species][chrom] = []
    
    # store list of CGIs that lift to that species
    inFile_lift = open(species+fileEnd_liftsToOtherSpecies, 'rt')
    for line in inFile_lift:
        splitLine = line.strip().split('\t')
        CGInameString = splitLine[3]
        chrom_in_hg38 = chromOfCGI_Dict[CGInameString]
        
        CGIs_lift_Dict[species][chrom_in_hg38].append(CGInameString)
        
    inFile_lift.close()
    
    # store list of CGIs that lift to that species AND don't overlap a feature
    inFile_lift_noFeature = open(species+fileEnd_liftsToOtherSpecies_andNoFeature+species+'.bed', 'rt')
    for line in inFile_lift_noFeature:
        splitLine = line.strip().split('\t')
        CGInameString = splitLine[3]
        chrom_in_hg38 = chromOfCGI_Dict[CGInameString]
        
        CGIs_lift_noFeature_Dict[species][chrom_in_hg38].append(CGInameString)
        
    inFile_lift.close()
    
print('checkpoint 1')

# re-open all CGIs in human
inAllCGIs_human = open('mergedCGIs_hg38.bed', 'rt')

# create blacklist to store CGIs that overlap a feature in another species
CGI_blacklist = {}
# CGI_blacklist[chrom] = [CGI1, CGI2, CGI3, ...]

for chrom in chromList:
    CGI_blacklist[chrom] = []

n = 0

for line in inAllCGIs_human:
    splitLine = line.strip().split('\t')
    chrom = splitLine[0]
    CGInameString = splitLine[3]

    for species in inSpeciesList:
        
        # if CGI overlaps a feature in another species, add it to the blacklist
        if CGInameString in CGIs_lift_Dict[species][chrom] and CGInameString not in CGIs_lift_noFeature_Dict[species][chrom]:
            if CGInameString not in CGI_blacklist[chrom]:
                CGI_blacklist[chrom].append(CGInameString)
                
    n += 1
    if n % 1000 == 0:
        print(n)
        
inAllCGIs_human.close()

print('checkpoint 2')

# read CGI list again and initialize each CGI in CGI_Dict if it's not blacklisted

CGI_Dict = {}
# CGI_Dict[chrom][CGInameString] = [status in species1, status in species2, ...]

for chrom in chromList:
    CGI_Dict[chrom] = {}
    
inAllCGIs_human = open('mergedCGIs_hg38.bed', 'rt')

n = 0

for line in inAllCGIs_human:
    splitLine = line.strip().split('\t')
    chrom = str(splitLine[0])
    CGInameString = splitLine[3]

    CGI_Dict[chrom][CGInameString] = [1,0,0,0,0,0,0,0,0]
            
    if 'hg38' in CGInameString:
        CGI_Dict[chrom][CGInameString][0] = 2
            
    for species in inSpeciesList:
        if CGInameString in CGIs_lift_noFeature_Dict[species][chrom]:
            CGI_Dict[chrom][CGInameString][SpeciesIndex_Dict[species]+1] = 1
        if species in CGInameString:
            CGI_Dict[chrom][CGInameString][SpeciesIndex_Dict[species]+1] = 2
                
    n += 1
    if n % 1000 == 0:
        print(n)
            
inAllCGIs_human.close()

print('checkpoint 3')

outFile = open('ak20240302_CGI_summary_acrossTree.txt', 'wt')

# print header
outFile.write('CGIname' +'\t'+ 'hg38' +'\t'+ '\t'.join(inSpeciesList)+'\t'+'phastCons'+'\n')

# write to output
for chrom in chromList:

    for CGInameString in CGI_Dict[chrom]:

        # convert to strings
        stringList = []
        for j in CGI_Dict[chrom][CGInameString]:
            stringList.append(str(j))
            
        # collect phastCons info
        if CGInameString in phastCons_Dict[chrom]:
            phastCons = 1
        else:
            phastCons = 0
    
        # print
        outFile.write(CGInameString +'\t'+ '\t'.join(stringList)+'\t'+str(phastCons)+'\n')
