# 4/4/22
# Modified heavily from 1/5/21 file "makeCategoryBedFiles_plusIntrons.py"
# This version uses NCBI gtf for biotype info, and uses it to combine NCBI_RefSeq/*.gtf and UCSC_RefSeq/*.gtf files
# (both from UCSC) into a set of 5 bed files with annotated regions
# Modified 4/19/22 to make a 7th output file for 5kb downstream of protein-coding promoters - removed on 7/24/22
# Modified 7/24/22 to consider chr sizes and not go beyond ends
# also ignore genes on patched regions since I won't incorporate them downstream

# Purpose: take 3 GTF files (described below), and output 6 bed files corresponding to the categories of genes from Bell & Vertino 2017:
# 1) protein-coding promoters (2kb upstream of TSS)
# 2) protein-coding exons
# 3) lncRNA exons + 2kb upstream
# 4) ncRNA exons + 2kb upstream
# 5) pseudogenes + 2kb upstream
# 6) unknown = anything without a biotype

# Usage: python makeCategoryBedFiles.py <species> <NCBI_RefSeq/*.gtf> <UCSC_RefSeq/*.gtf> <NCBI/*.gtf>
# <NCBI_RefSeq/*.gtf> = gtf file of NCBI sites but provided by UCSC
# <UCSC_RefSeq/*.gtf> = gtf file of NCBI transcripts remapped by UCSC
# <NCBI/*.gtf> = gtf provided by NCBI that contains "biotype" information

import sys

inSpecies = str(sys.argv[1])
inNCBI_GTF = open(sys.argv[2],'rt')
inUCSC_GTF = open(sys.argv[3],'rt')
inNCBI_biotype = open(sys.argv[4],'rt')
inChrSize = open(sys.argv[5],'rt')

outBed1 = open(inSpecies+"_1_proteinCoding_promoters_preMerge.bed",'wt')
outBed2 = open(inSpecies+"_2_proteinCoding_exons_preMerge.bed",'wt')
outBed3 = open(inSpecies+"_3_lncRNA_preMerge.bed",'wt')
outBed4 = open(inSpecies+"_4_ncRNA_preMerge.bed",'wt')
outBed5 = open(inSpecies+"_5_pseudogene_preMerge.bed",'wt')
outBed6 = open(inSpecies+"_6_unknown_preMerge.bed",'wt')

# I manually determined which of the 6 categories each feature type should belong to
# All 6 types get merged for the rest of the pipeline, but the assignment matters for determining whether to remove introns
FeatureType_Dict = {'antisense_RNA' : 4,'C_gene_segment' : 1,'C_region' : 1,'C_region_pseudogene' : 5,'D_gene_segment' : 1,'D_segment' : 1,'D_segment_pseudogene' : 5,'guide_RNA' : 4,'J_gene_segment' : 1,'J_segment' : 1,'J_segment_pseudogene' : 5,'lncRNA' : 3,'lnc_RNA' : 3,'miRNA' : 4,'misc_RNA' : 4,'mRNA' : 1,'ncRNA' : 4,'ncRNA_pseudogene' : 5,'other' : 1,'primary_transcript' : 1,'protein_coding' : 1,'pseudogene' : 5,'RNase_MRP_RNA' : 3,'RNase_P_RNA' : 4,'rRNA' : 4,'scRNA' : 4,'snoRNA' : 4,'snRNA' : 4,'SRP_RNA' : 4,'telomerase_RNA' : 4,'transcribed_pseudogene' : 5,'transcript' : 1,'tRNA' : 4,'vault_RNA' : 4,'V_gene_segment' : 1,'V_segment' : 1,'V_segment_pseudogene' : 5,'Y_RNA' : 4}

# make dictionary with chr lengths
ChrSize_Dict = {}
# ChrSize_Dict[chr] = length of chr

for line in inChrSize:
    splitLine = line.strip().split('\t')
    ChrSize_Dict[splitLine[0]] = int(splitLine[1])
inChrSize.close()

# Store biotype info from NCBI GTFs to sort genes into 6 categories

GeneName_Dict = {}
# GeneName_Dict[gene_id] = category (based on biotype translated with FeatureType_Dict)

for line in inNCBI_biotype:
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

                if splitField[0] == 'gene_biotype ' or splitField[0] == 'transcript_biotype ':
                    biotype = splitField[1]
                    category = FeatureType_Dict[biotype]

            # store in dictionary that converts gene name to category
            GeneName_Dict[geneID] = category

            geneID = ''
            biotype = ''
            category = ''
            
inNCBI_biotype.close()

#print(len(GeneName_Dict))
#for i in GeneName_Dict:
#    print(str(i)+'\t'+str(GeneName_Dict[i]))

# go through NCBI_RefSeq/*.gtf and print to appropriate output files using dictionary made from NCBI biotype GTF
for file in [inNCBI_GTF,inUCSC_GTF]:
    for line in file:
        splitLine = line.strip().split('\t')

        chrom = splitLine[0]
        if chrom in ChrSize_Dict:
        
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

                # calculate US 2kb
                if strand == '+':
                    US_start = start - 2000
                    US_end = start
                    if US_start < 0:
                        US_start = 0
                elif strand == '-':
                    US_start = end
                    US_end = end + 2000
                    if US_end > ChrSize_Dict[chrom]:
                        US_end = ChrSize_Dict[chrom]

                # only write to output if interval is greater than 0 bp
                # (fixes an error from bedtools merge making things go 1bp beyond the chr end)

                # protein-coding genes: write promoter to file 1 and 5kb downstream of promoter to file 7
                if category == 1 and US_end != 0 and US_start != US_end:
                    outBed1.write(str(chrom)+'\t'+str(US_start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')

                # same for categories 3, 4, 5
                if category == 3 and US_end != 0 and US_start != US_end:
                    outBed3.write(str(chrom)+'\t'+str(US_start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')
                if category == 4 and US_end != 0 and US_start != US_end:
                    outBed4.write(str(chrom)+'\t'+str(US_start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')
                if category == 5 and US_end != 0 and US_start != US_end:
                    outBed5.write(str(chrom)+'\t'+str(US_start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')

                # if no category, print entire record to file 6 (unknown)
                if category == 'NA':
                    if strand == '+' and US_start != end:
                        outBed6.write(str(chrom)+'\t'+str(US_start)+'\t'+str(end)+'\t'+str(geneID)+'\n')
                    if strand == '-' and start != US_end:
                        outBed6.write(str(chrom)+'\t'+str(start)+'\t'+str(US_end)+'\t'+str(geneID)+'\n')

            if splitLine[2] == 'exon':

                # protein-coding genes: write exon to file 2
                if category == 1 and start != end:
                    outBed2.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(geneID)+'\n')

                # same for categories 3,4,5
                if category == 3 and start != end:
                    outBed3.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(geneID)+'\n')
                if category == 4 and start != end:
                    outBed4.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(geneID)+'\n')
                if category == 5 and start != end:
                    outBed5.write(str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(geneID)+'\n')

            geneID = ''

inNCBI_GTF.close()
inUCSC_GTF.close()

outBed1.close()
outBed2.close()
outBed3.close()
outBed4.close()
outBed5.close()
outBed6.close()
