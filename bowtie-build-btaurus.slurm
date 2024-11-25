#!/bin/bash

#SBATCH --job-name=BOWTIE-BUILD-BOSTAURUS
#SBATCH --time=12:00:00
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --output=mpiOut-bowtie-build%j.out

bowtie-build /scratch/user/vignesh.bioenv/Novogene/Host_genomes/BosTaurus.fasta /scratch/user/vignesh.bioenv/Novogene/Host_genomes/Btaurus/genome_index
