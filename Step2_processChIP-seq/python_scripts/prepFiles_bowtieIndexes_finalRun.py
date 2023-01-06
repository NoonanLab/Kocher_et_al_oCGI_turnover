# 1/6/22
# Purpose: make files for making new bowtie indexes in 10 species for mapping Roller data
# Modified 7/11/22 to make indexes in new Gibbs storage space (/gpfs/gibbs/pi/noonan/ak2267/genomes) instead of scratch60

# python prepFiles_bowtieIndexes.py <genome>

import sys

genome = str(sys.argv[1])
outFile = open(genome+'_bowtieIndexes.batch','wt')

outFile.write("#!/bin/bash"+'\n'+"#SBATCH -p general"+'\n'+"#SBATCH -J "+genome+'_bowtieIndexes'+'\n'+"#SBATCH -n 1"+'\n'+"#SBATCH -c 10"+'\n'+"#SBATCH --mail-type=END"+'\n'+"#SBATCH --mail-user=acadia.kocher@yale.edu"+'\n'+'\n')
outFile.write("cd /gpfs/gibbs/pi/noonan/ak2267/genomes"+'\n'+"source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc"+'\n'+'\n')
outFile.write("module load Bowtie2/2.4.2-GCCcore-10.2.0" +'\n'+"module load SAMtools/1.12-GCCcore-10.2.0"+'\n'+'\n')

outFile.write("wget https://hgdownload.soe.ucsc.edu/goldenPath/"+genome+"/bigZips/"+genome+".fa.gz"+"\n")
outFile.write("bowtie2-build --threads 10 "+genome+".fa.gz "+genome)

# submit
# for i in $(ls *_bowtieIndexes.batch); do echo $i; sbatch $i; done
