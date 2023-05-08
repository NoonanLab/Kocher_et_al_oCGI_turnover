# 2/20/22
# Purpose: make files for making new bowtie indexes for mm39 and mm39_humanizedHs754

# python prepFiles_bowtieIndexes.py <genome>

import sys

genome = str(sys.argv[1])
outFile = open(genome+'_bowtieIndexes.batch','wt')

outFile.write("#!/bin/bash"+'\n'+"#SBATCH -p general"+'\n'+"#SBATCH -J "+genome+'_bowtieIndexes'+'\n'+"#SBATCH -n 1"+'\n'+"#SBATCH -c 10"+'\n'+"#SBATCH --mail-type=END"+'\n'+"#SBATCH --mail-user=acadia.kocher@yale.edu"+'\n'+'\n')
outFile.write("cd /home/ak2267/genomes"+'\n'+"source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc"+'\n'+'\n')
outFile.write("module load Bowtie2/2.4.2-GCCcore-10.2.0" +'\n'+"module load SAMtools/1.12-GCCcore-10.2.0"+'\n'+'\n')

outFile.write("bowtie2-build --threads 10 "+genome+".fa "+genome)

# submit
# for i in $(ls *_bowtieIndexes.batch); do echo $i; sbatch $i; done
