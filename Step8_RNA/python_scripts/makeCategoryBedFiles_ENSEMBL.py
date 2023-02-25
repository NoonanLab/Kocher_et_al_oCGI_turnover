# 5/19/22
# Modified heavily from 4/4/22 file "makeCategoryBedFiles_plusIntrons.py"

# Purpose: take ENSEMBL GTF file and output bed files corresponding to protein-coding promoters by GREAT rules (5kb up and 1kb down)

# Usage: python makeCategoryBedFiles_ENSEMBL.py <ENSEMBL GTF> <species in UCSC>

import sys

inENSEMBL_GTF = open(sys.argv[1],'rt')
inSpecies = str(sys.argv[2])

outBed = open(inSpecies+"_1_ENSEMBL_proteinCoding_promoters_preMerge.bed",'wt')

FeatureType_Dict = {'IG_C_gene' : 1,'IG_C_pseudogene' : 5,'IG_D_gene' : 1,'IG_D_pseudogene' : 5,'IG_J_gene' : 1,'IG_LV_gene' : 1,'IG_pseudogene' : 5,'IG_V_gene' : 1,'IG_V_pseudogene' : 5,'lncRNA' : 3,'miRNA' : 4,'misc_RNA' : 4,'Mt_rRNA' : 4,'Mt_tRNA' : 4,'polymorphic_pseudogene' : 5,'processed_pseudogene' : 5,'protein_coding' : 1,'pseudogene' : 5,'ribozyme' : 4,'rRNA' : 4,'scaRNA' : 4,'scRNA' : 4,'snoRNA' : 4,'snRNA' : 4,'sRNA' : 4,'TEC' : 6,'transcribed_processed_pseudogene' : 5,'transcribed_unitary_pseudogene' : 5,'transcribed_unprocessed_pseudogene' : 5,'translated_unprocessed_pseudogene' : 5,'TR_C_gene' : 1,'TR_D_gene' : 1,'TR_J_gene' : 1,'TR_J_pseudogene' : 5,'TR_V_gene' : 1,'TR_V_pseudogene' : 5,'unitary_pseudogene' : 5,'unprocessed_pseudogene' : 5,'vault_RNA' : 4,'Y_RNA' : 4}

# Store biotype info from ENSEMBL GTF to sort genes into categories

GeneName_Dict = {}
# GeneName_Dict[gene_id] = category (based on biotype translated with FeatureType_Dict)

for line in inENSEMBL_GTF:
    if line[0] != '#':
        #print(line)
        splitLine = line.strip().split('\t')

        if splitLine[2] == 'gene' or splitLine[2] == 'transcript':

            # get gene name and biotype for sorting into a category
            infoArray = splitLine[8].split('; ')
        
            # remove last ';' from array
            if infoArray[-1][-1] == ';':
                #print(infoArray)
                infoArray[-1] = infoArray[-1][:-1]
                #print(infoArray)

            for field in infoArray:
                #print(field)
                splitField = field.split('\"')
                #print(splitField)
                
                if splitField[0] == 'gene_id ':
                    geneID = splitField[1].split('_')[0]

                if splitField[0] == 'gene_biotype ':
                    biotype = splitField[1]
                    category = FeatureType_Dict[biotype]

            # store in dictionary that converts gene name to category
            GeneName_Dict[geneID] = category

            geneID = ''
            biotype = ''
            category = ''
            
inENSEMBL_GTF.close()

#print(len(GeneName_Dict))
#for i in GeneName_Dict:
#    print(str(i)+'\t'+str(GeneName_Dict[i]))

# go through GTF again and print promoters of protein-coding genes to output
inENSEMBL_GTF = open(sys.argv[1],'rt')

for line in inENSEMBL_GTF:
    if line[0] != '#':
        splitLine = line.strip().split('\t')
        chrom = splitLine[0]
        
        # GTFs are 1-based and I'm making 0-based bed files so subtract 1 from the start
        if int(splitLine[3]) > 0:
            start = int(splitLine[3]) - 1
        else:
            start = 0
        end = int(splitLine[4])
        strand = splitLine[6]

        # get gene name and biotype for sorting into a category
        infoArray = splitLine[8].split('; ')

        # remove last ';' from array
        if infoArray[-1][-1] == ';':
            #print(infoArray)
            infoArray[-1] = infoArray[-1][:-1]
            #print(infoArray)

        # go through fields to find gene_id
        for field in infoArray:
            #print(field)
            splitField = field.split('\"')
            #print(splitField)

            if splitField[0] == 'gene_id ':
                geneID = splitField[1]

        if geneID in GeneName_Dict:
            category = GeneName_Dict[geneID]
        else:
            category = 'NA'

        if splitLine[2] == 'transcript':    

            # calculate US 5kb & DS 1kb
            if strand == '+':
                US_start = start - 5000
                US_end = start + 1000
                if US_start < 0:
                    US_start = 0
            elif strand == '-':
                US_start = end - 1000
                US_end = end + 5000
                # MIGHT BE ISSUE HERE IF ANY GO OVER CHROM LENGTH

            # protein-coding genes: write promoter to file 1
            # also add 'chr' before chrom name
            #if category == 1 and US_end != 0 and :
            if category == 1 and len(chrom) < 3 and chrom != 'MT':
                outBed.write('chr'+str(chrom)+'\t'+str(US_start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')

        geneID = ''

inENSEMBL_GTF.close()

outBed.close()
