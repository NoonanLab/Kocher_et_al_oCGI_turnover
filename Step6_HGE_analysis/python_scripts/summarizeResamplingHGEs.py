# 2/15/23
# Purpose: summarize HGE resampling results

import sys

fileInfo = sys.argv[1]

inObserved = open(fileInfo + '_observed_HGEs_CGI_patterns.txt', 'rt')
inExpected = open(fileInfo + '_resampled_HGEs_CGI_patterns.txt', 'rt')

# STORE OBSERVED VALUES (i.e. counts from HGE set)

Observed_Dict = {}
# Observed_Dict['CGI_pattern'] = #

headerLine = 'yes'
for line in inObserved:
    splitLine = line.strip().split('\t')
    
    if headerLine == 'no':
        for i in range(0, len(splitLine)):
            Observed_Dict[header[i]] = int(splitLine[i])
            
    if headerLine == 'yes':
        header = splitLine
        headerLine = 'no'
        
inObserved.close()

#print(Observed_Dict)

# STORE EXPECTED VALUES (i.e. counts from resamplings)
    
Expected_Dict = {}
# Expected_Dict['CGI_pattern'] = [# in resampling 1, # in resampling 2, ...]

headerLine = 'yes'
for line in inExpected:
    splitLine = line.strip().split('\t')

    # if this line is not the header line, append to values array in Expected_Dict
    if headerLine == 'no':
    
        # append resampling value to array in Expected_Dict
        for i in range(0, len(splitLine)):
            resampleValue = 0
            if splitLine[i] != 'NA':
                resampleValue = int(splitLine[i])
            Expected_Dict[header[i]].append(resampleValue)
    
    # if this line is the header line, store as keys in Expected_Dict
    if headerLine == 'yes':
    
        header = splitLine
        headerLine = 'no'
        
        # initialize keys and value arrays in Expected_Dict
        for i in range(0, len(splitLine)):
            if header[i] not in Expected_Dict:
                Expected_Dict[header[i]] = []
    
#print(Expected_Dict)
    
inExpected.close()

#print(Expected_Dict['H'][0:9])
#print(len(Expected_Dict['H']))

# COMPARE OBSERVED AND EXPECTED AND WRITE TO OUTPUT
for CGI_pattern in Expected_Dict:
    
    observedValue = 0
    if CGI_pattern in Observed_Dict:
        observedValue = int(Observed_Dict[CGI_pattern])
    #print(len(Expected_Dict[CGI_pattern]))
    expectedValue = float(sum(Expected_Dict[CGI_pattern])) / float(len(Expected_Dict[CGI_pattern]))
    
    numberHigher = 0
    numberLower = 0
    
    for i in Expected_Dict[CGI_pattern]:
        if i > observedValue:
            numberHigher += 1
        if i < observedValue:
            numberLower += 1
            
    # calculate p value
    numberMoreExtreme = min(numberHigher, numberLower)
    if numberMoreExtreme == 0:
        numberMoreExtreme = 1
    
    p = float(numberMoreExtreme) / float(len(Expected_Dict[CGI_pattern]))
    
    # WRITE TO OUTPUT
    outArray = fileInfo.split('_') + [CGI_pattern, str(observedValue), str(expectedValue), str(p)]
    print('\t'.join(outArray))

