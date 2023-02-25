# 11/22/22
# Purpose: make table for R for each speciesPair x tissue x mark
# Updated 12/8/22 to incorporate 1) info on human orthology and 2) number of other peaks associated with each gene
# Updated 12/12/22 to work with regulatory domain files instead of bedtools closest

# FILES TO READ IN (hard-coded):
# 1) ENSEMBL 1:1 Orthologs: /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_mm39.txt 
# 2) linkage of oCGIs to regulatory domains in species A: rheMac10_mm39_brain_H3K4me3_speciesA_nearestUpstream.txt
# 3) linkage of oCGIs to regulatory domains in species B: rheMac10_mm39_brain_H3K4me3_speciesA_nearestDownstream.txt
# 4) read counts in species A - rep 1
# 5) read counts in species A - rep 2
# 6) read counts in species A - rep 3
# 7) read counts in species B - rep 1
# 8) read counts in species B - rep 2
# 9) read counts in species B - rep 3
# 10) human orthology file for speciesA
# 11) human orthology file for speciesB
# 12) intersection of THIS MARK peaks with regulatory domains in species A
# 13) intersection of THIS MARK peaks with regulatory domains in species B
# 14) intersection of MERGED intronic-intergenic peaks (all marks) with regulatory domains in species A
# 15) intersection of MERGED intronic-intergenic peaks (all marks) with regulatory domains in species B

# OUTPUT: table with the following columns:
# ENSG in A, ENSG in B
# Counts in A (rep 1, 2, 3), Counts in B (rep 1, 2, 3)
# Names of CGIs with Peaks in A, Names of CGIs with Peaks in B (these can be A-only OR B-only OR neither OR a mix)
# Whether CGIs are identical, Number of A-only CGIs with A-only peaks, Number of B-only CGIs with B-only peaks, CGI_summary
# Whether ENSG_A is 1:1 ortholog in human, Whether ENSG_B is 1:1 ortholog in human
# Number of same mark peaks in A, Number of same mark peaks in B, Number of all mark peaks in A, Number of all mark peaks in B

import sys

speciesPair = str(sys.argv[1])
speciesA = speciesPair.split('_')[0]
speciesB = speciesPair.split('_')[1]
tissue = str(sys.argv[2])
mark = str(sys.argv[3])
countFileEnd = str(sys.argv[4])

##### STEP 0: open input files

inOrthologs = open('/home/ak2267/genomes/ENSEMBL/1-1_Orthologs/'+ speciesPair +'.txt', 'rt')

inCGIs_RegulatoryDomains_A = open('overlappingRegulatoryDomains/'+ speciesPair +'_'+ tissue +'_'+ mark +'_speciesA.txt', 'rt')
inCGIs_RegulatoryDomains_B = open('overlappingRegulatoryDomains/'+ speciesPair +'_'+ tissue +'_'+ mark +'_speciesB.txt', 'rt')

inCounts_A1 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesA +'_'+ tissue +'_1'+countFileEnd)
inCounts_A2 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesA +'_'+ tissue +'_2'+countFileEnd)
inCounts_A3 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesA +'_'+ tissue +'_3'+countFileEnd)
inCounts_B1 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesB +'_'+ tissue +'_1'+countFileEnd)
inCounts_B2 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesB +'_'+ tissue +'_2'+countFileEnd)
inCounts_B3 = open('/home/ak2267/scratch60/Roller/bam_RNA/featureCounts/'+ speciesB +'_'+ tissue +'_3'+countFileEnd)

inHumanOrthology_A = open('/home/ak2267/genomes/ENSEMBL/1-1_Orthologs/hg38_'+speciesA+'.txt', 'rt')
inHumanOrthology_B = open('/home/ak2267/genomes/ENSEMBL/1-1_Orthologs/hg38_'+speciesB+'.txt', 'rt')

inThisMark_RegulatoryDomains_A = open('allPeaks/'+speciesA+'_'+tissue+'_'+mark+'_intersectRegulatoryDomains.txt', 'rt')
inThisMark_RegulatoryDomains_B = open('allPeaks/'+speciesB+'_'+tissue+'_'+mark+'_intersectRegulatoryDomains.txt', 'rt')

inAllMarks_RegulatoryDomains_A = open('allPeaks/'+speciesA+'_'+tissue+'_intersectRegulatoryDomains.txt', 'rt')
inAllMarks_RegulatoryDomains_B = open('allPeaks/'+speciesB+'_'+tissue+'_intersectRegulatoryDomains.txt', 'rt')

##### STEP 1: store dictionaries of orthologs

# first between speciesA and speciesB

Ortholog_Dict = {}
Ortholog_Dict['A'] = {}
# Ortholog_Dict['A'][ENSG_A] = [ENSG_B_1, ENSG_B_2, ...]
Ortholog_Dict['B'] = {}
# Ortholog_Dict['B'][ENSG_B] = [ENSG_A_1, ENSG_A_2, ...]

for line in inOrthologs:
    if line[0] != 'G':
        splitLine = line.strip().split('\t')
        ENSG_A = splitLine[0]
        ENSG_B = splitLine[4]
        
        # store in A dict
        if ENSG_A not in Ortholog_Dict['A']:
            Ortholog_Dict['A'][ENSG_A] = []
        if ENSG_B not in Ortholog_Dict['A'][ENSG_A]:
            Ortholog_Dict['A'][ENSG_A].append(ENSG_B)
        
        # store in B dict
        if ENSG_B not in Ortholog_Dict['B']:
            Ortholog_Dict['B'][ENSG_B] = []
        if ENSG_A not in Ortholog_Dict['B'][ENSG_B]:
            Ortholog_Dict['B'][ENSG_B].append(ENSG_A)
inOrthologs.close()

# then between human and speciesA, and human and speciesB
HumanOrtholog_Dict = {}
HumanOrtholog_Dict['humanToA'] = {}
# HumanOrtholog_Dict['humanToA'][ENSG_human] = [ENSG_A_1, ENSG_A_2, ...]
HumanOrtholog_Dict['humanToB'] = {}
# HumanOrtholog_Dict['humanToB'][ENSG_human] = [ENSG_B_1, ENSG_B_2, ...]
HumanOrtholog_Dict['AtoHuman'] = {}
# HumanOrtholog_Dict['AtoHuman'][ENSG_A] = [ENSG_human_1, ENSG_human_2, ...]
HumanOrtholog_Dict['BtoHuman'] = {}
# HumanOrtholog_Dict['BtoHuman'][ENSG_B] = [ENSG_human_1, ENSG_human_2, ...]

for line in inHumanOrthology_A:
    if line[0] != 'G':
        splitLine = line.strip().split('\t')
        ENSG_human = splitLine[0]
        ENSG_A = splitLine[4]
        
        # store in HumanOrtholog_Dict (both ways)
        if ENSG_human not in HumanOrtholog_Dict['humanToA']:
            HumanOrtholog_Dict['humanToA'][ENSG_human] = []
        if ENSG_A not in HumanOrtholog_Dict['humanToA'][ENSG_human]:
            HumanOrtholog_Dict['humanToA'][ENSG_human].append(ENSG_A)
            
        # store in HumanOrtholog_Dict (both ways)
        if ENSG_A not in HumanOrtholog_Dict['AtoHuman']:
            HumanOrtholog_Dict['AtoHuman'][ENSG_A] = []
        if ENSG_human not in HumanOrtholog_Dict['AtoHuman'][ENSG_A]:
            HumanOrtholog_Dict['AtoHuman'][ENSG_A].append(ENSG_human)
inHumanOrthology_A.close()

for line in inHumanOrthology_B:
    if line[0] != 'G':
        splitLine = line.strip().split('\t')
        ENSG_human = splitLine[0]
        ENSG_B = splitLine[4]
        
        # store in HumanOrtholog_Dict (both ways)
        if ENSG_human not in HumanOrtholog_Dict['humanToB']:
            HumanOrtholog_Dict['humanToB'][ENSG_human] = []
        if ENSG_B not in HumanOrtholog_Dict['humanToB'][ENSG_human]:
            HumanOrtholog_Dict['humanToB'][ENSG_human].append(ENSG_B)
            
        # store in HumanOrtholog_Dict (both ways)
        if ENSG_B not in HumanOrtholog_Dict['BtoHuman']:
            HumanOrtholog_Dict['BtoHuman'][ENSG_B] = []
        if ENSG_human not in HumanOrtholog_Dict['BtoHuman'][ENSG_B]:
            HumanOrtholog_Dict['BtoHuman'][ENSG_B].append(ENSG_human)
inHumanOrthology_B.close()


#for ENSG_A in HumanOrtholog_Dict['AtoHuman']:
 #   print(ENSG_A)
  #  print(HumanOrtholog_Dict['AtoHuman'][ENSG_A])
   # print(len(HumanOrtholog_Dict['AtoHuman'][ENSG_A]))
    #if len(HumanOrtholog_Dict['AtoHuman'][ENSG_A]) == 1:
     #   ENSG_human = HumanOrtholog_Dict['AtoHuman'][ENSG_A][0]
      #  print(ENSG_human)
       # print(HumanOrtholog_Dict['humanToA'][ENSG_human])
        #print(len(HumanOrtholog_Dict['humanToA'][ENSG_human]))
    

##### STEP 2: link orthologs to species-specific oCGIs with species-specific peaks

# inCGIs_RegulatoryDomains_A = open('overlappingRegulatoryDomains/'+ speciesPair +'_'+ tissue +'_'+ mark +'_speciesA.txt', 'rt')
# inCGIs_RegulatoryDomains_B = open('overlappingRegulatoryDomains/'+ speciesPair +'_'+ tissue +'_'+ mark +'_speciesB.txt', 'rt')

Nearest_Dict = {}
Nearest_Dict['A'] = {}
# Nearest_Dict['A'][ENSG_A] = [CGIname1, CGIname2, ...]
Nearest_Dict['B'] = {}
# Nearest_Dict['B'][ENSG_B] = [CGIname1, CGIname2, ...]

FileList = [inCGIs_RegulatoryDomains_A, inCGIs_RegulatoryDomains_B]
SpeciesList = ['A', 'B']

for i in range(0,2):
    for line in FileList[i]:
        splitLine = line.strip().split('\t')
        CGIname = splitLine[3]
        ENSG = splitLine[7]
        
        if ENSG not in Nearest_Dict[SpeciesList[i]]:
            Nearest_Dict[SpeciesList[i]][ENSG] = []
        if CGIname not in Nearest_Dict[SpeciesList[i]][ENSG]:
            Nearest_Dict[SpeciesList[i]][ENSG].append(CGIname)
                
##### STEP 2.5: link orthologs to their intronic/intergenic peaks (all marks)

#inThisMark_RegulatoryDomains_A = open('allPeaks/'+speciesA+'_'+tissue+'_'+mark+'_intersectRegulatoryDomains.txt', 'rt')
#inThisMark_RegulatoryDomains_B = open('allPeaks/'+speciesB+'_'+tissue+'_'+mark+'intersectRegulatoryDomains.txt', 'rt')

#inAllMarks_RegulatoryDomains_A = open('allPeaks/'+speciesA+'_'+tissue+'_intersectRegulatoryDomains.txt', 'rt')
#inAllMarks_RegulatoryDomains_B = open('allPeaks/'+speciesB+'_'+tissue+'_intersectRegulatoryDomains.txt', 'rt')

FileList = [inThisMark_RegulatoryDomains_A, inThisMark_RegulatoryDomains_B, inAllMarks_RegulatoryDomains_A, inAllMarks_RegulatoryDomains_B]
IndexList = ['ThisMark_A', 'ThisMark_B', 'AllMarks_A', 'AllMarks_B']

NearestAllPeaks_Dict = {}
NearestAllPeaks_Dict['ThisMark_A'] = {}
# NearestAllPeaks_Dict['ThisMark_A'][ENSG_A] = [peakName1, peakName2, ...]
NearestAllPeaks_Dict['ThisMark_B'] = {}
# NearestAllPeaks_Dict['ThisMark_B'][ENSG_B] = [peakName1, peakName2, ...]
NearestAllPeaks_Dict['AllMarks_A'] = {}
# NearestAllPeaks_Dict['AllMarks_A'][ENSG_A] = [peakName1, peakName2, ...]
NearestAllPeaks_Dict['AllMarks_B'] = {}
# NearestAllPeaks_Dict['AllMarks_B'][ENSG_B] = [peakName1, peakName2, ...]

for i in range(0,4):
    #print(i)
    for line in FileList[i]:
        splitLine = line.strip().split('\t')
        peakName = splitLine[3]
        ENSG = splitLine[7]
        if ENSG != '.':
            if ENSG not in NearestAllPeaks_Dict[IndexList[i]]:
                NearestAllPeaks_Dict[IndexList[i]][ENSG] = []
            if peakName not in NearestAllPeaks_Dict[IndexList[i]][ENSG]:
                NearestAllPeaks_Dict[IndexList[i]][ENSG].append(peakName)
            
#for i in NearestAllPeaks_Dict['ThisMark_A']:
 #   print(str(i) +'\t'+ str(NearestAllPeaks_Dict['ThisMark_A'][i]))

##### STEP 3: store counts for each ortholog in both species

if '.counts' in countFileEnd:
    countIndex = 6
elif 'ReadsPerGene.' in countFileEnd:
    countIndex = 3

Counts_Dict = {}
Counts_Dict['A'] = {}
# Counts_Dict['A'][ENSG_A] = [counts in rep1, counts in rep2, counts in rep3]
Counts_Dict['B'] = {}
# Counts_Dict['B'][ENSG_B] = [counts in rep1, counts in rep2, counts in rep3]

FileList = [inCounts_A1, inCounts_A2, inCounts_A3, inCounts_B1, inCounts_B2, inCounts_B3]
SpeciesList = ['A', 'A', 'A', 'B', 'B', 'B']
RepList = [1, 2, 3, 1, 2, 3]

for i in range(0, 6):
    for line in FileList[i]:
        if line[0] != 'N' and line[0] != '#' and line[0] != 'G':
            splitLine = line.strip().split('\t')
            ENSG = splitLine[0]
            counts = int(splitLine[countIndex])
            
            if ENSG not in Counts_Dict[SpeciesList[i]]:
                Counts_Dict[SpeciesList[i]][ENSG] = [0, 0, 0]
            Counts_Dict[SpeciesList[i]][ENSG][RepList[i] - 1] = counts
    FileList[i].close()

##### STEP 4: write to output - only for 1:1 orthologs (skip 1:many, many:1, and many:many)

# ENSG in A, ENSG in B
# Counts in A (rep 1, 2, 3), Counts in B (rep 1, 2, 3)
# Names of CGIs with Peaks in A, Names of CGIs with Peaks in B (these can be A-only OR B-only OR neither OR a mix)
# Whether CGIs are identical, Number of A-only CGIs with A-only peaks, Number of B-only CGIs with B-only peaks, CGI_summary
header = []
header += ['ENSG_A', 'ENSG_B']
header += ['HumanOrtholog_A', 'HumanOrtholog_B']
header += ['Counts_A1', 'Counts_A2', 'Counts_A3', 'Counts_B1', 'Counts_B2', 'Counts_B3']
header += ['CGIs_A', 'CGIs_B', 'SameCGIs', 'Num_A_only', 'Num_B_only', 'CGI_summary']
header += ['Peaks_ThisMark_A', 'Peaks_ThisMark_B', 'Peaks_AllMarks_A', 'Peaks_AllMarks_B']
print('\t'.join(header))

# Ortholog_Dict['A'][ENSG_A] = [ENSG_B_1, ENSG_B_2, ...]
# Ortholog_Dict['B'][ENSG_B] = [ENSG_A_1, ENSG_A_2, ...]

# HumanOrtholog_Dict['humanToA'][ENSG_human] = [ENSG_A_1, ENSG_A_2, ...]
# HumanOrtholog_Dict['humanToB'][ENSG_human] = [ENSG_B_1, ENSG_B_2, ...]
# HumanOrtholog_Dict['AtoHuman'][ENSG_A] = [ENSG_human_1, ENSG_human_2, ...]
# HumanOrtholog_Dict['BtoHuman'][ENSG_B] = [ENSG_human_1, ENSG_human_2, ...]

# Nearest_Dict['A'][ENSG_A] = [CGIname1, CGIname2, ...]
# Nearest_Dict['B'][ENSG_B] = [CGIname1, CGIname2, ...]

# Counts_Dict['A'][ENSG_A] = [counts in rep1, counts in rep2, counts in rep3]
# Counts_Dict['B'][ENSG_B] = [counts in rep1, counts in rep2, counts in rep3]

# NearestAllPeaks_Dict['ThisMark_A'][ENSG_A] = [peakName1, peakName2, ...]
# NearestAllPeaks_Dict['ThisMark_B'][ENSG_B] = [peakName1, peakName2, ...]
# NearestAllPeaks_Dict['AllMarks_A'][ENSG_A] = [peakName1, peakName2, ...]
# NearestAllPeaks_Dict['AllMarks_B'][ENSG_B] = [peakName1, peakName2, ...]

countOrthologNotInGTF = 0

for ENSG_A in Ortholog_Dict['A']:

    # test if ENSG_A goes with only one ENSG_B
    if len(Ortholog_Dict['A'][ENSG_A]) == 1:
        ENSG_B = Ortholog_Dict['A'][ENSG_A][0]
        # test if ENSG_B goes with only one ENSG_A
        if len(Ortholog_Dict['B'][ENSG_B]) == 1:
            
            # if so, start collecting info for writing to output
            # ENSG_A, ENSG_B
            outputArray = []
            outputArray.append(ENSG_A)
            outputArray.append(ENSG_B)
            
            # see if gene is 1:1 ortholog with human gene for both A and B
            # HumanOrtholog_A, HumanOrtholog_B
            HumanOrtholog_A = 0
            HumanOrtholog_B = 0
            
            if ENSG_A in HumanOrtholog_Dict['AtoHuman']:
                if len(HumanOrtholog_Dict['AtoHuman'][ENSG_A]) == 1:
                    ENSG_human = HumanOrtholog_Dict['AtoHuman'][ENSG_A][0]
                    # test if 1:1 in other direction (human to A)
                    if len(HumanOrtholog_Dict['humanToA'][ENSG_human]) == 1:
                        HumanOrtholog_A = 1
            if ENSG_B in HumanOrtholog_Dict['BtoHuman']:
                if len(HumanOrtholog_Dict['BtoHuman'][ENSG_B]) == 1:
                    ENSG_human = HumanOrtholog_Dict['BtoHuman'][ENSG_B][0]
                    # test if 1:1 in other direction (human to B)
                    if len(HumanOrtholog_Dict['humanToB'][ENSG_human]) == 1:
                        HumanOrtholog_B = 1
            outputArray.append(HumanOrtholog_A)
            outputArray.append(HumanOrtholog_B)
            
            #print(str(HumanOrtholog_A) +'\t'+ str(HumanOrtholog_B))
            
            # collect counts
            # Counts_A1, Counts_A2, Counts_A3, Counts_B1, Counts_B2, Counts_B3
            if ENSG_A in Counts_Dict['A']:
                outputArray += Counts_Dict['A'][ENSG_A]
            else:
                outputArray += ['NA', 'NA', 'NA']
                countOrthologNotInGTF += 1
            if ENSG_B in Counts_Dict['B']:
                outputArray += Counts_Dict['B'][ENSG_B]
            else:
                outputArray += ['NA', 'NA', 'NA']
                
            # collect CGInames
            # CGIs_A, CGIs_B'
            if ENSG_A in Nearest_Dict['A']:
                outputArray.append(';'.join(Nearest_Dict['A'][ENSG_A]))
            else:
                outputArray.append('NA')
            if ENSG_B in Nearest_Dict['B']:
                outputArray.append(';'.join(Nearest_Dict['B'][ENSG_B]))
            else:
                outputArray.append('NA')
            
            # make summary info
            # SameCGIs, Num_A_only, Num_B_only, CGI_summary
            sameCGIs = '0'
            Num_A_only = 'NA'
            Num_B_only = 'NA'
            CGI_summary = 'NA'
            
            if ENSG_A in Nearest_Dict['A'] and ENSG_B in Nearest_Dict['B']:
                if sorted(Nearest_Dict['A'][ENSG_A]) == sorted(Nearest_Dict['B'][ENSG_B]):
                    sameCGIs = '1'
                    Num_A_only = 0
                    Num_B_only = 0
                    for i in Nearest_Dict['A'][ENSG_A]:
                        if 'a_' in i:
                            Num_A_only += 1
                        if 'b_' in i:
                            Num_B_only += 1
                    if Num_A_only > 0 and Num_B_only == 0:
                        CGI_summary = 'A'
                    if Num_A_only == 0 and Num_B_only > 0:
                        CGI_summary = 'B'
                    if Num_A_only > 0 and Num_B_only > 0:
                        CGI_summary = 'mix'
             
            outputArray.append(sameCGIs)
            outputArray.append(Num_A_only)
            outputArray.append(Num_B_only)            
            outputArray.append(CGI_summary)
            
            # info on other genes in the region
            # Peaks_ThisMark_A, Peaks_ThisMark_B, Peaks_AllMarks_A, Peaks_AllMarks_B
            Peaks_ThisMark_A = 0
            Peaks_ThisMark_B = 0
            Peaks_AllMarks_A = 0
            Peaks_AllMarks_B = 0
                
            if ENSG_A in NearestAllPeaks_Dict['ThisMark_A']:
                Peaks_ThisMark_A = len(NearestAllPeaks_Dict['ThisMark_A'][ENSG_A])
                Peaks_AllMarks_A = len(NearestAllPeaks_Dict['AllMarks_A'][ENSG_A])
            if ENSG_B in NearestAllPeaks_Dict['ThisMark_B']:
                Peaks_ThisMark_B = len(NearestAllPeaks_Dict['ThisMark_B'][ENSG_B])
                Peaks_AllMarks_B = len(NearestAllPeaks_Dict['AllMarks_B'][ENSG_B])
                
            outputArray.append(Peaks_ThisMark_A)
            outputArray.append(Peaks_ThisMark_B)
            outputArray.append(Peaks_AllMarks_A)
            outputArray.append(Peaks_AllMarks_B)
            
            outputArrayAsStrings = []
            for i in outputArray:
                outputArrayAsStrings.append(str(i))
            
            print('\t'.join(outputArrayAsStrings))
    
#print(countOrthologNotInGTF)

