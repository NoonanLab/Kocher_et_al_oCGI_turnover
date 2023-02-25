# 2/12/23
# Purpose: take output of RNA_resamplingTest.R and make summary across all species pairs x tissues x marks

import sys

inSpeciesPair = sys.argv[1]
inTissue = sys.argv[2]
inMark = sys.argv[3]

fileStem = inSpeciesPair + '_' + inTissue + '_' + inMark

# read and store observed values
inObserved = open(fileStem + '_observedMedians.txt')

for line in inObserved:
    splitLine = line.strip().split()
    A_obs = float(splitLine[0])
    B_obs = float(splitLine[1])
inObserved.close()

# read and store expected values
inExpected_A = open(fileStem + '_A_resamplingMedians.txt')
inExpected_B = open(fileStem + '_B_resamplingMedians.txt')

A_exp_array = []
B_exp_array = []

for line in inExpected_A:
    exp = float(line.strip())
    A_exp_array.append(exp)
inExpected_A.close()

for line in inExpected_B:
    exp = float(line.strip())
    B_exp_array.append(exp)
inExpected_B.close()

# calculate mean of resampling medians
A_exp = sum(A_exp_array) / len(A_exp_array)
B_exp = sum(B_exp_array) / len(B_exp_array)

# calculate number of resampling medians that are more extreme than observed
# only test A being higher and B being lower

A_higher_count = 0
for i in A_exp_array:
    if i > A_obs:
        A_higher_count += 1

B_lower_count = 0
for i in B_exp_array:
    if i < B_obs:
        B_lower_count += 1
        
# replace with 1 if none were more extreme
if A_higher_count == 0:
    A_higher_count = 1
if B_lower_count == 0:
    B_lower_count = 1
        
# calculate p values
p_A = A_higher_count / len(A_exp_array)
p_B = B_lower_count / len(B_exp_array)   

# write to output
outArray_A = [inSpeciesPair, inTissue, inMark, 'A', str(A_obs), str(A_exp), str(p_A)]
print('\t'.join(outArray_A))

outArray_B = [inSpeciesPair, inTissue, inMark, 'B', str(B_obs), str(B_exp), str(p_B)]
print('\t'.join(outArray_B)) 
