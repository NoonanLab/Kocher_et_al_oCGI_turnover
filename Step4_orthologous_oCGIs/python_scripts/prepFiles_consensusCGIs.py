# 7/30/22
# Purpose: create jobFile for generating consensus CGIs for:
# 1) all 28 two-way combinations of Roller species
# 2) all 3 two-way combinations of Noonan species

speciesArrayRoller = ['rheMac10','calJac4','mm39','rn7','susScr11','canFam6','felCat9','equCab3']
speciesArrayNoonan = ['hg19','rheMac2','mm9']

SpeciesCombos = []

# store filePaths and fileEnds
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs'
pathToCGIs = '/home/ak2267/genomes/CGI/UCSC_AL/'
fileEndCGIs = '_CGIsAL.bed'
FilesToExclude_Dict = {}
FilesToExclude_Dict['filePath'] = {'RefSeq':'/home/ak2267/genomes/RefSeq/featureAnnotations/','FANTOM':'/home/ak2267/genomes/FANTOM/FANTOM_TSS_','blacklist':'/home/ak2267/genomes/blacklist/blacklist_'}
FilesToExclude_Dict['fileEnd'] = {'RefSeq':'_allFeatures.bed','FANTOM':'.bed','blacklist':'.bed'}

# make list of species combinations (pairwise)
for j in range(0,len(speciesArrayRoller)):
    for k in range(j+1,len(speciesArrayRoller)):
        speciesString = str(speciesArrayRoller[j]) +'_'+ str(speciesArrayRoller[k])
        SpeciesCombos.append(speciesString)
for j in range(0,len(speciesArrayNoonan)):
    for k in range(j+1,len(speciesArrayNoonan)):
        speciesString = str(speciesArrayNoonan[j]) +'_'+ str(speciesArrayNoonan[k])
        SpeciesCombos.append(speciesString)

for speciesCombo in SpeciesCombos:
    #print(speciesCombo)
    speciesA = speciesCombo.split('_')[0]
    speciesB = speciesCombo.split('_')[1]

    outputString = []

    #########
    # STEP 0: append info on where to work and source files
    outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
    outputString.append('mkdir '+speciesCombo+'; cd '+speciesCombo)

    #########
    # STEP 1: make files with 4 columns and CGI names (as "peak" names), with 'a' or 'b' preceding in each species
    #outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' '+pathToCGIs+speciesA+fileEndCGIs+' > speciesA_CGI_4col.bed')
    #outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' '+pathToCGIs+speciesB+fileEndCGIs+' > speciesB_CGI_4col.bed')

    #########
    # ALT STEP 1: first merge CGIs within 200bp, then make files with 4 columns and CGI names (as "peak" names), with 'a' or 'b' preceding in each species
    outputString.append('bedtools merge -d 200 -i '+pathToCGIs+speciesA+fileEndCGIs+' | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' > speciesA_CGI_4col.bed')
    outputString.append('bedtools merge -d 200 -i '+pathToCGIs+speciesB+fileEndCGIs+' | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' > speciesB_CGI_4col.bed')

    #########
    # STEP 2: get 'orphan' CpG islands in each species that liftover to the other species

    speciesWithFANTOMandBlacklist = ['hg19','mm39','mm9']
    
    # store file paths for restricting to non-feature (aka orphan) CGIs in both species
    if speciesA in speciesWithFANTOMandBlacklist:
        filesToExcludeA = [FilesToExclude_Dict['filePath']['RefSeq']+speciesA+FilesToExclude_Dict['fileEnd']['RefSeq'],
                          FilesToExclude_Dict['filePath']['FANTOM']+speciesA+FilesToExclude_Dict['fileEnd']['FANTOM'],
                          FilesToExclude_Dict['filePath']['blacklist']+speciesA+FilesToExclude_Dict['fileEnd']['blacklist']]
    else:
        filesToExcludeA = [FilesToExclude_Dict['filePath']['RefSeq']+speciesA+FilesToExclude_Dict['fileEnd']['RefSeq']]
    if speciesB in speciesWithFANTOMandBlacklist:
        filesToExcludeB = [FilesToExclude_Dict['filePath']['RefSeq']+speciesB+FilesToExclude_Dict['fileEnd']['RefSeq'],
                          FilesToExclude_Dict['filePath']['FANTOM']+speciesB+FilesToExclude_Dict['fileEnd']['FANTOM'],
                          FilesToExclude_Dict['filePath']['blacklist']+speciesB+FilesToExclude_Dict['fileEnd']['blacklist']]
    else:
        filesToExcludeB = [FilesToExclude_Dict['filePath']['RefSeq']+speciesB+FilesToExclude_Dict['fileEnd']['RefSeq']]
        
    hg38FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg38'+FilesToExclude_Dict['fileEnd']['RefSeq'],
                          FilesToExclude_Dict['filePath']['FANTOM']+'hg38'+FilesToExclude_Dict['fileEnd']['FANTOM'],
                          FilesToExclude_Dict['filePath']['blacklist']+'hg38'+FilesToExclude_Dict['fileEnd']['blacklist']]
    hg19FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg19'+FilesToExclude_Dict['fileEnd']['RefSeq'],
                          FilesToExclude_Dict['filePath']['FANTOM']+'hg19'+FilesToExclude_Dict['fileEnd']['FANTOM'],
                          FilesToExclude_Dict['filePath']['blacklist']+'hg19'+FilesToExclude_Dict['fileEnd']['blacklist']]


    # lift both species to human

    # remove features in both species
    outputString.append('bedtools intersect -v -a speciesA_CGI_4col.bed -b '+' -b '.join(filesToExcludeA)+' > speciesA_CGI_4col_noFeature.bed')
    outputString.append('bedtools intersect -v -a speciesB_CGI_4col.bed -b '+' -b '.join(filesToExcludeB)+' > speciesB_CGI_4col_noFeature.bed')
    
    # ROLLER SPECIES

    if speciesA in speciesArrayRoller and speciesB in speciesArrayRoller:
        chainAtoHg38 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg38.over.chain.gz'
        chainBtoHg38 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg38.over.chain.gz'
        chainHg38toA = '/home/ak2267/genomes/chain/hg38To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
        chainHg38toB = '/home/ak2267/genomes/chain/hg38To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

        # lift to human
        outputString.append('liftOver -minMatch=0.3 speciesA_CGI_4col_noFeature.bed '+chainAtoHg38+' speciesA_CGI_4col_noFeature_LOhg38.bed unMapped')
        outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature.bed '+chainBtoHg38+' speciesB_CGI_4col_noFeature_LOhg38.bed unMapped')

        # lift back and restrict list in human to only those that lift back to the same CGI
        outputString.append('liftOver -minMatch=0.3 speciesA_CGI_4col_noFeature_LOhg38.bed '+chainHg38toA+' speciesA_CGI_4col_noFeature_liftedBackToA.bed unMapped')
        outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature_LOhg38.bed '+chainHg38toB+' speciesB_CGI_4col_noFeature_liftedBackToB.bed unMapped')
        outputString.append('bedtools intersect -wao -a speciesA_CGI_4col_noFeature_liftedBackToA.bed -b speciesA_CGI_4col_noFeature.bed > speciesA_checkMappingBack.txt')
        outputString.append('bedtools intersect -wao -a speciesB_CGI_4col_noFeature_liftedBackToB.bed -b speciesB_CGI_4col_noFeature.bed > speciesB_checkMappingBack.txt')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBack.txt > speciesA_CGI_4col_noFeature_mapBack_inAcoord.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBack.txt > speciesB_CGI_4col_noFeature_mapBack_inBcoord.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_CGI_4col_noFeature_LOhg38.bed speciesA_CGI_4col_noFeature_mapBack_inAcoord.bed > speciesA_CGI_4col_noFeature_LOhg38_mapBack.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_CGI_4col_noFeature_LOhg38.bed speciesB_CGI_4col_noFeature_mapBack_inBcoord.bed > speciesB_CGI_4col_noFeature_LOhg38_mapBack.bed')
        
        # merge, remove features in human, sort file names to be a_# then b_#
        outputString.append('cat speciesA_CGI_4col_noFeature_LOhg38_mapBack.bed speciesB_CGI_4col_noFeature_LOhg38_mapBack.bed > cat.bed')
        outputString.append('sort -k1,1 -k2,2n cat.bed > sort.bed')
        outputString.append('bedtools merge -c 4 -o collapse -i sort.bed > merge_hg38.bed')
        outputString.append('bedtools intersect -v -a merge_hg38.bed -b '+' -b '.join(hg38FilesToExclude)+' > merge_hg38_noFeature.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/sortRegionNames.py merge_hg38_noFeature.bed > merge_hg38_noFeature_sortedNames.bed')

        # lift back out to speciesA and speciesB
        outputString.append('liftOver -minMatch=0.3 merge_hg38_noFeature_sortedNames.bed '+chainHg38toA+' mergedCGIs_liftToA.bed unMapped')
        outputString.append('liftOver -minMatch=0.3 merge_hg38_noFeature_sortedNames.bed '+chainHg38toB+' mergedCGIs_liftToB.bed unMapped')

        # lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
        outputString.append('liftOver -minMatch=0.3 mergedCGIs_liftToA.bed '+chainAtoHg38+' mergedCGIs_liftToA_backToHg38.bed unMapped')
        outputString.append('liftOver -minMatch=0.3 mergedCGIs_liftToB.bed '+chainBtoHg38+' mergedCGIs_liftToB_backToHg38.bed unMapped')
        outputString.append('bedtools intersect -wao -a mergedCGIs_liftToA_backToHg38.bed -b merge_hg38_noFeature_sortedNames.bed > speciesA_checkMappingBackInHg38.txt')
        outputString.append('bedtools intersect -wao -a mergedCGIs_liftToB_backToHg38.bed -b merge_hg38_noFeature_sortedNames.bed > speciesB_checkMappingBackInHg38.txt')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBackInHg38.txt > speciesA_merged_thatMapBackToHg38_inHg38.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBackInHg38.txt > speciesB_merged_thatMapBackToHg38_inHg38.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToA.bed speciesA_merged_thatMapBackToHg38_inHg38.bed > merged_speciesA_mapBack.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToB.bed speciesB_merged_thatMapBackToHg38_inHg38.bed > merged_speciesB_mapBack.bed')

        # remove features
        outputString.append('bedtools intersect -v -a merged_speciesA_mapBack.bed -b '+' -b '.join(filesToExcludeA)+' > mergedCGIs_liftToA_noFeature.bed')
        outputString.append('bedtools intersect -v -a merged_speciesB_mapBack.bed -b '+' -b '.join(filesToExcludeB)+' > mergedCGIs_liftToB_noFeature.bed')

        #########
        # STEP 3:
        # restrict all three lists (in A, in B, and in hg38) to those present in both A and B
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToA_noFeature.bed mergedCGIs_liftToB_noFeature.bed > '+speciesCombo+'_speciesA_reconciledCGI.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToB_noFeature.bed mergedCGIs_liftToA_noFeature.bed > '+speciesCombo+'_speciesB_reconciledCGI.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py merge_hg38_noFeature_sortedNames.bed mergedCGIs_liftToA_noFeature.bed > intermediateHg38_reconciledFile.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intermediateHg38_reconciledFile.bed mergedCGIs_liftToB_noFeature.bed > '+speciesCombo+'_hg38_reconciledCGI.bed')
        
    # NOONAN SPECIES
    elif speciesA in speciesArrayNoonan and speciesB in speciesArrayNoonan:
        chainAtoHg19 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg19.over.chain.gz'
        chainBtoHg19 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg19.over.chain.gz'
        chainHg19toA = '/home/ak2267/genomes/chain/hg19To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
        chainHg19toB = '/home/ak2267/genomes/chain/hg19To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

        # lift to human
        if speciesA == 'hg19':
            outputString.append('cp speciesA_CGI_4col_noFeature.bed speciesA_CGI_4col_noFeature_LOhg19.bed')
            outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature.bed '+chainBtoHg19+' speciesB_CGI_4col_noFeature_LOhg19.bed unMapped')
        elif speciesA != 'hg19':
            outputString.append('liftOver -minMatch=0.3 speciesA_CGI_4col_noFeature.bed '+chainAtoHg19+' speciesA_CGI_4col_noFeature_LOhg19.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature.bed '+chainBtoHg19+' speciesB_CGI_4col_noFeature_LOhg19.bed unMapped')

        # lift back and restrict list in human to only those that lift back to the same CGI
        if speciesA == 'hg19':
            outputString.append('cp speciesA_CGI_4col_noFeature_LOhg19.bed speciesA_CGI_4col_noFeature_liftedBackToA.bed')
            outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature_LOhg19.bed '+chainHg19toB+' speciesB_CGI_4col_noFeature_liftedBackToB.bed unMapped')
        elif speciesA != 'hg19':
            outputString.append('liftOver -minMatch=0.3 speciesA_CGI_4col_noFeature_LOhg19.bed '+chainHg19toA+' speciesA_CGI_4col_noFeature_liftedBackToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_CGI_4col_noFeature_LOhg19.bed '+chainHg19toB+' speciesB_CGI_4col_noFeature_liftedBackToB.bed unMapped')
        outputString.append('bedtools intersect -wao -a speciesA_CGI_4col_noFeature_liftedBackToA.bed -b speciesA_CGI_4col_noFeature.bed > speciesA_checkMappingBack.txt')
        outputString.append('bedtools intersect -wao -a speciesB_CGI_4col_noFeature_liftedBackToB.bed -b speciesB_CGI_4col_noFeature.bed > speciesB_checkMappingBack.txt')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBack.txt > speciesA_CGI_4col_noFeature_mapBack_inAcoord.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBack.txt > speciesB_CGI_4col_noFeature_mapBack_inBcoord.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_CGI_4col_noFeature_LOhg19.bed speciesA_CGI_4col_noFeature_mapBack_inAcoord.bed > speciesA_CGI_4col_noFeature_LOhg19_mapBack.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_CGI_4col_noFeature_LOhg19.bed speciesB_CGI_4col_noFeature_mapBack_inBcoord.bed > speciesB_CGI_4col_noFeature_LOhg19_mapBack.bed')
        
        # merge, remove features in human, sort file names to be a_# then b_#
        outputString.append('cat speciesA_CGI_4col_noFeature_LOhg19_mapBack.bed speciesB_CGI_4col_noFeature_LOhg19_mapBack.bed > cat.bed')
        outputString.append('sort -k1,1 -k2,2n cat.bed > sort.bed')
        outputString.append('bedtools merge -c 4 -o collapse -i sort.bed > merge_hg19.bed')
        outputString.append('bedtools intersect -v -a merge_hg19.bed -b '+' -b '.join(hg19FilesToExclude)+' > merge_hg19_noFeature.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/sortRegionNames.py merge_hg19_noFeature.bed > merge_hg19_noFeature_sortedNames.bed')

        # lift back out to speciesA and speciesB
        if speciesA == 'hg19':
            outputString.append('cp merge_hg19_noFeature_sortedNames.bed mergedCGIs_liftToA.bed')
            outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toB+' mergedCGIs_liftToB.bed unMapped')            
        elif speciesA != 'hg19':
            outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toA+' mergedCGIs_liftToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toB+' mergedCGIs_liftToB.bed unMapped')

       # lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
        if speciesA == 'hg19':
            outputString.append('cp mergedCGIs_liftToA.bed mergedCGIs_liftToA_backToHg19.bed')
            outputString.append('liftOver -minMatch=0.3 mergedCGIs_liftToB.bed '+chainBtoHg19+' mergedCGIs_liftToB_backToHg19.bed unMapped')            
        elif speciesA != 'hg19':
            outputString.append('liftOver -minMatch=0.3 mergedCGIs_liftToA.bed '+chainAtoHg19+' mergedCGIs_liftToA_backToHg19.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 mergedCGIs_liftToB.bed '+chainBtoHg19+' mergedCGIs_liftToB_backToHg19.bed unMapped')
        outputString.append('bedtools intersect -wao -a mergedCGIs_liftToA_backToHg19.bed -b merge_hg19_noFeature_sortedNames.bed > speciesA_checkMappingBackInHg19.txt')
        outputString.append('bedtools intersect -wao -a mergedCGIs_liftToB_backToHg19.bed -b merge_hg19_noFeature_sortedNames.bed > speciesB_checkMappingBackInHg19.txt')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBackInHg19.txt > speciesA_merged_thatMapBackToHg19_inHg19.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBackInHg19.txt > speciesB_merged_thatMapBackToHg19_inHg19.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToA.bed speciesA_merged_thatMapBackToHg19_inHg19.bed > merged_speciesA_mapBack.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToB.bed speciesB_merged_thatMapBackToHg19_inHg19.bed > merged_speciesB_mapBack.bed')

        # remove features
        outputString.append('bedtools intersect -v -a merged_speciesA_mapBack.bed -b '+' -b '.join(filesToExcludeA)+' > mergedCGIs_liftToA_noFeature.bed')
        outputString.append('bedtools intersect -v -a merged_speciesB_mapBack.bed -b '+' -b '.join(filesToExcludeB)+' > mergedCGIs_liftToB_noFeature.bed')


        #########
        # STEP 3:
        # restrict all three lists (in A, in B, and in hg19) to those present in both A and B
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToA_noFeature.bed mergedCGIs_liftToB_noFeature.bed > '+speciesCombo+'_speciesA_reconciledCGI.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_liftToB_noFeature.bed mergedCGIs_liftToA_noFeature.bed > '+speciesCombo+'_speciesB_reconciledCGI.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py merge_hg19_noFeature_sortedNames.bed mergedCGIs_liftToA_noFeature.bed > intermediateHg19_reconciledFile.bed')
        outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intermediateHg19_reconciledFile.bed mergedCGIs_liftToB_noFeature.bed > '+speciesCombo+'_hg19_reconciledCGI.bed')
        
    print(' ; '.join(outputString))
