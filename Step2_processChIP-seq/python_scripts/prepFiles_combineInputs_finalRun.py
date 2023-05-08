# 1/7/22
# Purpose: make files for batch submission for combining input files
# since several come from ArrayExpress separated into two files

# cd /home/ak2267/project/Roller/ChIP/batch
# python prepFiles/prepFiles_combineInputs.py

fileArray = ['canFam6_muscle_input_3','equCab3_testis_input_2','felCat9_muscle_input_2','rheMac10_muscle_input_2','rheMac10_muscle_input_3','susScr11_muscle_input_1']

for fileString in fileArray:

    outFile = open(str(fileString)+'.batch','wt')

    outFile.write("#!/bin/bash"+'\n'+"#SBATCH -p general"+'\n'+"#SBATCH -J "+fileString+'\n'+"#SBATCH -n 1"+'\n'+"#SBATCH -c 10"+'\n'+"#SBATCH --mail-type=END"+'\n'+"#SBATCH --mail-user=acadia.kocher@yale.edu"+'\n'+'\n')
    outFile.write("cd /gpfs/gibbs/pi/noonan/ak2267/Roller/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc"+'\n')
    outFile.write("module load SAMtools/1.12-GCCcore-10.2.0"+'\n')
    outFile.write("samtools merge "+str(fileString)+"_bowtie2_filtered.bam "+str(fileString)+"a_bowtie2_filtered.bam "+str(fileString)+"b_bowtie2_filtered.bam"+'\n')
    outFile.close()
    
# submit
# for i in $(ls *.batch); do echo $i; sbatch $i; done
