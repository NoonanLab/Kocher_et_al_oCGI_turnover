# 3/19/2024
# Purpose: take output of faCount and intersections with feature annotations
# make summary tables for use in R
# modified from 3/5/2024 file to include information on N bases in the sequence

import sys

inFaCount = open(sys.argv[1], 'rt')
inNoFeatures = open(sys.argv[2], 'rt')
inRmsk = open(sys.argv[3], 'rt')

outTable = open(sys.argv[4], 'wt')

# store info on feature overlap in oCGI_List
oCGI_List = []

for line in inNoFeatures:
    splitLine = line.strip().split('\t')
    name = splitLine[0] +':'+ splitLine[1] +'-'+ splitLine[2]
    oCGI_List.append(name)
inNoFeatures.close()

# store info on repeat overlap from intersection with rmsk files
Rmsk_Dict = {}
# Rmsk_Dict[name] = # bases overlapping a repeat

for line in inRmsk:
    splitLine = line.strip().split('\t')
    name = splitLine[0] +':'+ splitLine[1] +'-'+ splitLine[2]
    if name not in Rmsk_Dict:
        Rmsk_Dict[name] = 0
    bp = int(splitLine[8])
    Rmsk_Dict[name] += bp
inRmsk.close()

# write header to output file
outTable.write('Interval' +'\t'+ 'Model_based' +'\t'+ 'UCSC_AL' +'\t'+ 'UCSC_Standard' +'\t'+ 'Length' +'\t'+ 'Rmsk_bp' +'\t'+ 'N_number' +'\t'+ 'GC_percent' +'\t'+ 'CpG_number' +'\t'+ 'Obs_Exp_Ratio' +'\t'+ 'oCGI' + '\n')

# open faCount file, manipulate info, and print to output table
n = 0
for line in inFaCount:
    if line[0] != '#' and line[0] != 't':
        splitLine = line.strip().split('\t')
        
        name = splitLine[0].split('::')[1]
        whichList = splitLine[0].split('::')[0]
        
        # make CGI_specificity into 3 columns
        CGI_specificity = [0,0,0]
        if 'Model_based' in whichList:
            CGI_specificity[0] = 1
        if 'UCSC_AL' in whichList:
            CGI_specificity[1] = 1
        if 'UCSC_Standard' in whichList:
            CGI_specificity[2] = 1
        
        # collect info
        length = int(splitLine[1])
        GC_pct = (int(splitLine[3]) + int(splitLine[4])) / int(splitLine[1]) * 100
        CpG_num = int(splitLine[7])
        N = int(splitLine[6])
        Rmsk_bp = Rmsk_Dict[name]
        
        #if int(splitLine[3]) == 0 or int(splitLine[4]) == 0:
         #   print(line)
         
        # print count if site has 0 Cs or 0 Gs (very small number of sites)
        if int(splitLine[3]) == 0 or int(splitLine[4]) == 0:
            n += 1
            #print(n)
            
            oE = 'NA'
            
        else:
        
            # observed / expected CpG ratio
            # CpG_num * length / (C_num * G_num)
            oE = CpG_num * length / (int(splitLine[3]) * int(splitLine[4]))
            
        # oCGI status
        oCGI_status = 0
        if name in oCGI_List:
            oCGI_status = 1
            
        # convert output into strings in an array
        outputArray = [str(name), str(CGI_specificity[0]), str(CGI_specificity[1]), str(CGI_specificity[2]), str(length), str(Rmsk_bp), str(N), str(GC_pct), str(CpG_num), str(oE), str(oCGI_status)]
        
        # write to output
        outTable.write('\t'.join(outputArray) + '\n')