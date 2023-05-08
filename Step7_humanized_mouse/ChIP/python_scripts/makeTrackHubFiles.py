# 10/29/22
# Purpose: make trackDb.txt for hosting hs754 ChIP data

#multiplexGroups = ['A', 'B', 'C', 'D', 'E']
multiplexGroups = ['D', 'A', 'E', 'B', 'C']

Mark_Dict = {'A':['H3K4me3', 'CTCF', 'H3K27me3'], 'B':['H3K4me3'], 'C':['H3K27me3'], 'D':['H3K27ac'], 'E':['H3K27ac']}
# Mark_Dict[multiplexGroup] = [mark1, mark2, ...]

PeakType_Dict = {'H3K27ac':'narrow', 'H3K4me3':'narrow', 'H3K27me3':'broad', 'CTCF':'narrow'}
# PeakType_Dict[mark] = broad/narrow

Rep_Dict = {'A':3, 'B':2, 'C':2, 'D':3, 'E':2}
# Rep_Dict[multiplexGroup] = # reps

Time_Dict = {'A':'e17.5', 'B':'e11.5', 'C':'e11.5', 'D':'e17.5', 'E':'e11.5'}
# Time_Dict[multiplexGroup] = timePoint

Color_Dict = {'H3K27ac':'59,179,0', 'H3K4me3':'255,153,51', 'H3K27me3':'89, 134, 207', 'CTCF':'22, 79, 171', 'input':'209,209,209'} 
# Color_Dict[mark] = RGB code for color

MaxView_Dict = {'H3K27ac':'1', 'H3K4me3':'2', 'H3K27me3':'5', 'CTCF':'1', 'input':'1'} 
# MaxView_Dict[mark] = max value for y axis in visualization

priorityCount = 0

# print entry for file marking humanized region (mouse replaced region)
print('track humanizedRegion')
print('type bigBed')
print('shortLabel humanizedRegion')
print('longLabel humanizedRegion')
print('visibility full')
print('color 35, 35, 35')
print('bigDataUrl mouseCoord_mm39.bb')
print('priority '+str(priorityCount)+'\n')

priorityCount += 1

# print entry for file marking humanized region (showing full humanized region extending 390bp)
print('track humanizedRegion_390bp')
print('type bigBed')
print('shortLabel humanizedRegion_390bp')
print('longLabel humanizedRegion_390bp')
print('visibility full')
print('color 35, 35, 35')
print('bigDataUrl humanCoord_mm39.bb')
print('priority '+str(priorityCount)+'\n')

priorityCount += 1

for multiplexGroup in multiplexGroups:

    repNum = Rep_Dict[multiplexGroup]
    timePoint = Time_Dict[multiplexGroup]

    for mark in Mark_Dict[multiplexGroup]:
        RGBcolor = Color_Dict[mark]
        peakType = PeakType_Dict[mark]
        maxView = MaxView_Dict[mark]
        
        for genotype in ['WT', 'HUM']:   
        
            for rep in range(1, repNum + 1):
                replicate = str(rep)
        
                print('track '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'_signal')
                print('type bigWig 0 '+str(maxView))
                print('shortLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
                print('longLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
                print('visibility full')
                print('color '+RGBcolor)
                print('windowingFunction mean')
                print('smoothingWindow 5')
                print('bigDataUrl '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'.bw')
                print('priority '+str(priorityCount)+'\n')

                priorityCount += 1
                
                print('track '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'_peaks')
                print('type bigBed')
                print('shortLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
                print('longLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
                print('visibility dense')
                print('color '+RGBcolor)
                print('bigDataUrl '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'_'+peakType+'Peak.bb')
                print('priority '+str(priorityCount)+'\n')
                
                priorityCount += 1
                
    mark = 'input'
    RGBcolor = Color_Dict[mark]
    maxView = MaxView_Dict[mark]
    
    for genotype in ['WT', 'HUM']:
    
        for rep in range(1, repNum + 1):
            replicate = str(rep)
        
            print('track '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'_signal')
            print('type bigWig 0 1')
            print('shortLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
            print('longLabel '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate)
            print('visibility full')
            print('color '+RGBcolor)
            print('windowingFunction mean')
            print('smoothingWindow 5')
            print('bigDataUrl '+multiplexGroup+'_'+timePoint+'_'+mark+'_'+genotype+'_'+replicate+'.bw')
            print('priority '+str(priorityCount)+'\n')

            priorityCount += 1

            


