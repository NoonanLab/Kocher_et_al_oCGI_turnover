# 9/21/22
# Purpose: incorporate all data generated in 220912_CGI_activityHeatmap.sh into single tables for 
# hg38 plus the 8 Roller species, incorporating Noonan data when it exists

# python makeSingleSpeciesSummaryTables_CGIcentric.py <species>

Roller_speciesArray = ['rheMac10','calJac4','mm39','rn7','susScr11','canFam6','felCat9','equCab3']
Noonan_speciesArray = ['hg38','rheMac10','mm39']

basePath = '/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/'

import sys
species = str(sys.argv[1])

##### STEP 0: get intronic-intergenic CGIs
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/${species}_CGIsAL_noFeatures.bed
##### STEP 1: intersect with ChIP-seq peaks
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPeaks/${species}_Roller_${tissue}_${mark}_intersection.txt
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPeaks/${species}_Noonan_brainAc_${timePoint}_intersection.txt
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPeaks/${species}_Noonan_brainMe2_${timePoint}_intersection.txt
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPeaks/${species}_Noonan_limbAc_${timePoint}_intersection.txt
##### STEP 2: run featureCounts to count reads in all peaks
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts/${species}_${tissue}_${mark}.signal				 # Roller
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts/${species}_${tissue}_${mark}_${timePoint}.signal   # Noonan
##### STEP 3: count CpGs in CGIs and peaks using faCount
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_CGIsAL.faCount			    			# CGIs
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_${tissue}_${mark}.faCount			    # Roller
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_${tissue}_${mark}_${timePoint}.faCount   # Noonan
##### STEP 4: intersect phastCons with CGIs and peaks - tricky because phastCons is in hg38
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/${species}_CGIsAL_intersectPhastConsInHg38.txt
##### STEP 5: intersect phastBias with CGIs and peaks - tricky because phastBias is in hg38
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastBias/${species}_CGIsAL_intersectPhastBiasInHg38.txt


CGI_Dict = {}
# CGI_Dict[CGIname] = [length, CpG num, phastCons max LOD, sum of LOD, sum of bp with LOD, 
# number of phastBias tracts, bp overlapping phastBias tract, total length of phastBias tract, coord in Roller genomes, liftsToHg38]

CGI_ActivityDict = {}
# CGI_ActivityDict[CGIname][tissue_mark (Roller) / tissueMark_timePoint (Noonan)] = [peak1, peak2, etc]
# just initialize outer dictionary here - inner will be later

WhetherCGIinNoonan_Dict = {}
# WhetherCGIinNoonan_Dict[CGIname] = 1/0

inCGI = open(basePath+species+'_CGIsAL_noFeatures.bed','rt')
for line in inCGI:
    splitLine = line.strip().split('\t')
    CGIname = splitLine[3]
    if CGIname not in CGI_Dict:
        CGI_Dict[CGIname] = [0,0,0,0,0,0,0,0,0,0]
        # CGI_Dict[CGIname] = [length, CpG num, phastCons max LOD, sum of LOD, sum of bp with LOD, 
        # number of phastBias tracts, bp overlapping phastBias tract, total length of phastBias tract, coord in Roller genomes, liftsToHg38]
        CGI_ActivityDict[CGIname] = {} # initialize for later
    CGIlength = int(splitLine[2]) - int(splitLine[1])
    CGI_Dict[CGIname][0] += CGIlength
    CGI_Dict[CGIname][8] = splitLine[0]+':'+splitLine[1]+'-'+splitLine[2]
    WhetherCGIinNoonan_Dict[CGIname] = 0
inCGI.close()

inCGI_faCount = open(basePath+'faCount/'+species+'_CGIsAL.faCount','rt')
for line in inCGI_faCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        CGIname = splitLine[0].split('::')[0]
        CpGnum = int(splitLine[7])
        CGI_Dict[CGIname][1] += CpGnum
inCGI_faCount.close()

inCGI_phastCons = open(basePath+'intersectPhastCons/'+species+'_CGIsAL_intersectPhastConsInHg38.txt','rt')
for line in inCGI_phastCons:
    splitLine = line.strip().split('\t')
    CGIname = splitLine[3]
    # if CGIname is in this file, it means it lifts to hg38 regardless of overlaps
    CGI_Dict[CGIname][9] = 1
    if splitLine[4] != '.':
        LOD = int(splitLine[7].split('=')[1])
        overlapBp = int(splitLine[8])
        if LOD > CGI_Dict[CGIname][2]:
           CGI_Dict[CGIname][2] = LOD
        CGI_Dict[CGIname][3] += LOD
        CGI_Dict[CGIname][4] += overlapBp
inCGI_phastCons.close()

inCGI_phastBias = open(basePath+'intersectPhastBias/'+species+'_CGIsAL_intersectPhastBiasInHg38.txt','rt')
for line in inCGI_phastBias:
    splitLine = line.strip().split('\t')
    CGIname = splitLine[3]
    if splitLine[4] != '.':
        overlap = int(splitLine[8])
        phastBiasLength = int(splitLine[6]) - int(splitLine[5])
        CGI_Dict[CGIname][5] += 1
        CGI_Dict[CGIname][6] += overlap
        CGI_Dict[CGIname][7] += phastBiasLength 
        #### PROBLEM - some CGIs overlap more than one phastBias tract - see felCat9
        # --> output the number so that you can filter later if you want
        # may want to also filter for the opposite, which is the number of CGIs per phastBias tract
        # leave for later!

# ROLLER DATA
Roller_tissueArray = ['brain','liver','muscle','testis']
Roller_markArray = ['H3K27ac','H3K4me3','H3K4me1']

Peak_Dict = {}
# Peak_Dict[peakName] = [length, CpGnum, RPM, RPKM]

# Add to CGI_ActivityDict
# CGI_ActivityDict[CGIname][tissue_mark (Roller) / tissueMark_timePoint (Noonan)] = [peak1, peak2, etc]

if species in Roller_speciesArray:
    for tissue in Roller_tissueArray:
        for mark in Roller_markArray:
        
            inPeakIntersection = open(basePath+'intersectPeaks/'+species+'_Roller_'+tissue+'_'+mark+'_intersection.txt','rt')
            for line in inPeakIntersection:
                splitLine = line.strip().split('\t')
                CGIname = splitLine[3]
                if CGIname not in CGI_ActivityDict:
                    CGI_ActivityDict[CGIname] = {}
                if tissue+'_'+mark not in CGI_ActivityDict[CGIname]:
                    CGI_ActivityDict[CGIname][tissue+'_'+mark] = []
                
                peakName = splitLine[7]
                if peakName != '.':
                    peakLength = int(splitLine[6]) - int(splitLine[5])
                
                    # FOUND CGIS CAN HAVE MULTIPLE PEAKS - output column saying how many
                    CGI_ActivityDict[CGIname][tissue+'_'+mark].append(peakName)
                    
                    if peakName not in Peak_Dict:
                        Peak_Dict[peakName] = [0,0,0,0]
                    Peak_Dict[peakName][0] = peakLength                           
            inPeakIntersection.close()
         
            # faCount
            inFaCount = open(basePath+'faCount/'+species+'_'+tissue+'_'+mark+'.faCount')
            for line in inFaCount:
                if line[0] != '#':
                    splitLine = line.strip().split('\t')
                    peakName = splitLine[0].split('::')[0]
                    CpGnum = int(splitLine[7])
                    if peakName not in Peak_Dict:
                        Peak_Dict[peakName] = [0,0,0,0]
                    Peak_Dict[peakName][1] = CpGnum
            inFaCount.close()
            
            # featureCounts
            inFeatureCounts = open(basePath+'featureCounts/'+species+'_'+tissue+'_'+mark+'.signal')
            for line in inFeatureCounts:
                splitLine = line.strip().split('\t')
                peakName = splitLine[0]
                RPM = splitLine[1]
                RPKM = splitLine[2]
                if peakName not in Peak_Dict:
                        Peak_Dict[peakName] = [0,0,0,0]
                Peak_Dict[peakName][2] = RPM
                Peak_Dict[peakName][3] = RPKM
            inFeatureCounts.close()
            
# NOONAN DATA
Noonan_tissueMarkArray = ['brainAc','brainMe2','limbAc']
Noonan_timePointArray = [0,1,2,3]

Noonan_SpeciesTranslator = {'hg38':'hg19','rheMac10':'rheMac2','mm39':'mm9'}
Noonan_MarkTranslator = {'brainAc':'ac','brainMe2':'me2','limbAc':'ac'}
Noonan_TissueTranslator = {'brainAc':'brain','brainMe2':'brain','limbAc':'limb'}
Noonan_TimePointTranslator = {}
Noonan_TimePointTranslator['hg19'] = {}
Noonan_TimePointTranslator['rheMac2'] = {}
Noonan_TimePointTranslator['mm9'] = {}
Noonan_TimePointTranslator['hg19']['brainAc'] = {0:'CS16', 1:'CS23', 2:'F2F', 3:'F2O'}
Noonan_TimePointTranslator['rheMac2']['brainAc'] = {0:'NA', 1:'e55', 2:'e79F', 3:'e79O'}
Noonan_TimePointTranslator['mm9']['brainAc'] = {0:'e11', 1:'e14', 2:'17F', 3:'17O'}
Noonan_TimePointTranslator['hg19']['brainMe2'] = {0:'CS16', 1:'CS23', 2:'F2F', 3:'F2O'}
Noonan_TimePointTranslator['rheMac2']['brainMe2'] = {0:'NA', 1:'e55', 2:'e79F', 3:'e79O'}
Noonan_TimePointTranslator['mm9']['brainMe2'] = {0:'e11', 1:'e14', 2:'17F', 3:'17O'}
Noonan_TimePointTranslator['hg19']['limbAc'] = {0:'E33', 1:'E41', 2:'E44', 3:'E47'}
Noonan_TimePointTranslator['rheMac2']['limbAc'] = {0:'NA', 1:'e31', 2:'e36', 3:'NA'}
Noonan_TimePointTranslator['mm9']['limbAc'] = {0:'e10.5', 1:'e11.5', 2:'e12.5', 3:'e13.5'}


if species in Noonan_speciesArray:
    for tissueMark in Noonan_tissueMarkArray:
        for timePoint in Noonan_timePointArray:
        
            # skip rhesus tissues with no data
            if Noonan_TimePointTranslator[Noonan_SpeciesTranslator[species]][tissueMark][timePoint] != 'NA':
        
                inPeakIntersection = open(basePath+'intersectPeaks/'+Noonan_SpeciesTranslator[species]+'_Noonan_'+tissueMark+'_'+Noonan_TimePointTranslator[Noonan_SpeciesTranslator[species]][tissueMark][timePoint]+'_intersection.txt','rt')
                for line in inPeakIntersection:
                    splitLine = line.strip().split('\t')
                    CGIname = splitLine[3]
                    WhetherCGIinNoonan_Dict[CGIname] = 1
                    
                    if tissueMark+'_'+str(timePoint) not in CGI_ActivityDict[CGIname]:
                        CGI_ActivityDict[CGIname][tissueMark+'_'+str(timePoint)] = []
                
                    peakName = tissueMark+'_'+str(timePoint)+'_'+splitLine[7]
                    if splitLine[7] != '.':
                        peakLength = int(splitLine[6]) - int(splitLine[5])
                
                        # FOUND CGIS CAN HAVE MULTIPLE PEAKS - output column saying how many
                        CGI_ActivityDict[CGIname][tissueMark+'_'+str(timePoint)].append(peakName)
                    
                        if peakName not in Peak_Dict:
                            Peak_Dict[peakName] = [0,0,0,0]
                        Peak_Dict[peakName][0] = peakLength                           
                inPeakIntersection.close()
         
                # faCount
                inFaCount = open(basePath+'faCount/'+Noonan_SpeciesTranslator[species]+'_'+Noonan_TissueTranslator[tissueMark]+'_'+Noonan_MarkTranslator[tissueMark]+'_'+Noonan_TimePointTranslator[Noonan_SpeciesTranslator[species]][tissueMark][timePoint]+'.faCount')
                for line in inFaCount:
                    if line[0] != '#':
                        splitLine = line.strip().split('\t')
                        peakName = tissueMark+'_'+str(timePoint)+'_'+splitLine[0].split('::')[0]
                        CpGnum = int(splitLine[7])
                        if peakName not in Peak_Dict:
                            Peak_Dict[peakName] = [0,0,0,0]
                        Peak_Dict[peakName][1] = CpGnum
                inFaCount.close()
            
                # featureCounts
                inFeatureCounts = open(basePath+'featureCounts/'+Noonan_SpeciesTranslator[species]+'_'+Noonan_TissueTranslator[tissueMark]+'_'+Noonan_MarkTranslator[tissueMark]+'_'+str(timePoint)+'.signal')
                for line in inFeatureCounts:
                    splitLine = line.strip().split('\t')
                    peakName = tissueMark+'_'+str(timePoint)+'_'+splitLine[0]
                    RPM = splitLine[1]
                    RPKM = splitLine[2]
                    if peakName not in Peak_Dict:
                            Peak_Dict[peakName] = [0,0,0,0]
                    Peak_Dict[peakName][2] = RPM
                    Peak_Dict[peakName][3] = RPKM
                inFeatureCounts.close()
            
# prepare headers
CGI_header = ['CGI_name','CGI_coord','CGI_length','CGI_CpGnum','CGI_liftsToHg38','CGI_phastConsMax','CGI_phastConsSum','CGI_phastConsBp','CGI_phastBiasNumber','CGI_phastBiasBpOverlap','CGI_phastBiasBpTotal']
Roller_header = []
for tissue in Roller_tissueArray:
    for mark in Roller_markArray:
        for headerEntry in ['peak','peakNumber','peakLength','CpGnum','avgRPM','avgRPKM']:
            Roller_header.append(headerEntry+'.'+tissue+'.'+mark)

# go through species and print to output
outFile = open(basePath+'summaryTables_singleSpeciesCGIcentric/'+species+'_summaryTable_singleSpeciesCGIcentric.txt','wt')
    
# collect output info on CGI
header = CGI_header
if species in Roller_speciesArray:
    header += Roller_header
if species in Noonan_speciesArray:
    # generate Noonan header, leaving out rhesus combinations with no data
    Noonan_header = ['LiftsToNoonanGenomes']
    for tissueMark in Noonan_tissueMarkArray:
        for timePoint in Noonan_timePointArray:
            if (species == 'rheMac10' and (tissueMark == 'brainAc' or tissueMark == 'brainMe2') and timePoint == 0) or (species == 'rheMac10' and tissueMark == 'limbAc' and (timePoint == 0 or timePoint == 3)):
                continue
            else:
                for headerEntry in ['peak','peakNumber','peakLength','CpGnum','avgRPM','avgRPKM']:
                    Noonan_header.append(headerEntry+'.'+tissueMark+'.'+str(timePoint))
    header += Noonan_header
outFile.write('\t'.join(header)+'\n')

# collect info for each CGI

for CGIname in CGI_Dict:
    outputArray = []
    outputArray.append(CGIname)
    outputArray.append(CGI_Dict[CGIname][8])
    outputArray.append(CGI_Dict[CGIname][0])
    outputArray.append(CGI_Dict[CGIname][1])
    outputArray.append(CGI_Dict[CGIname][9])
    outputArray.append(CGI_Dict[CGIname][2])
    outputArray.append(CGI_Dict[CGIname][3])
    outputArray.append(CGI_Dict[CGIname][4])
    outputArray.append(CGI_Dict[CGIname][5])
    outputArray.append(CGI_Dict[CGIname][6])
    outputArray.append(CGI_Dict[CGIname][7])

    # output Roller data
    if species in Roller_speciesArray:
        for tissue in Roller_tissueArray:
            for mark in Roller_markArray:
                if len(CGI_ActivityDict[CGIname][tissue+'_'+mark]) > 0:
                    outputArray.append(1)
                    peakNumber = len(CGI_ActivityDict[CGIname][tissue+'_'+mark])
                    outputArray.append(peakNumber)
                    # calculate across peaks
                    peakLength = 0
                    peakCpGnum = 0
                    peakTotalRPM = 0
                    peakTotalRPKM = 0
                    for peakName in CGI_ActivityDict[CGIname][tissue+'_'+mark]:
                        peakLength += Peak_Dict[peakName][0]
                        peakCpGnum += Peak_Dict[peakName][1]
                        peakTotalRPM += float(Peak_Dict[peakName][2])
                        peakTotalRPKM += float(Peak_Dict[peakName][3])
                    outputArray.append(peakLength)
                    outputArray.append(peakCpGnum)
                    outputArray.append(peakTotalRPM / peakNumber)
                    outputArray.append(peakTotalRPKM / peakNumber)
                else:
                    outputArray += [0] * 6
                    
    # output Noonan data (if it is a Noonan species)
    if species in Noonan_speciesArray:
        if WhetherCGIinNoonan_Dict[CGIname] == 1:
            outputArray.append(1)
        else:
            outputArray.append(0)
            
        for tissueMark in Noonan_tissueMarkArray:
            for timePoint in Noonan_timePointArray:
                if (species == 'rheMac10' and (tissueMark == 'brainAc' or tissueMark == 'brainMe2') and timePoint == 0) or (species == 'rheMac10' and tissueMark == 'limbAc' and (timePoint == 0 or timePoint == 3)):
                    continue
                else:
                    if CGIname in CGI_ActivityDict:
                        if tissueMark+'_'+str(timePoint) in CGI_ActivityDict[CGIname]:
                            if len(CGI_ActivityDict[CGIname][tissueMark+'_'+str(timePoint)]) > 0:
                                outputArray.append(1)
                                peakNumber = len(CGI_ActivityDict[CGIname][tissueMark+'_'+str(timePoint)])
                                outputArray.append(peakNumber)
                                # calculate across peaks
                                peakLength = 0
                                peakCpGnum = 0
                                peakTotalRPM = 0
                                peakTotalRPKM = 0
                                for peakName in CGI_ActivityDict[CGIname][tissueMark+'_'+str(timePoint)]:
                                    peakLength += Peak_Dict[peakName][0]
                                    peakCpGnum += Peak_Dict[peakName][1]
                                    peakTotalRPM += float(Peak_Dict[peakName][2])
                                    peakTotalRPKM += float(Peak_Dict[peakName][3])
                                outputArray.append(peakLength)
                                outputArray.append(peakCpGnum)
                                outputArray.append(peakTotalRPM / peakNumber)
                                outputArray.append(peakTotalRPKM / peakNumber)
                            else:
                                outputArray += [0] * 6
                        else:
                            outputArray += [0] * 6
                    else:
                        outputArray += [0] * 6


    outputAsStrings = []
    for i in outputArray:
        outputAsStrings.append(str(i))
        
    outFile.write('\t'.join(outputAsStrings)+'\n')
outFile.close()

# CGI_ActivityDict[CGIname][tissue_mark (Roller) / tissueMark_timePoint (Noonan)] = [peak1, peak2, etc]


