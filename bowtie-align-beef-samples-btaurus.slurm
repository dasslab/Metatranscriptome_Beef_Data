#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=FPKM-job           # job name
#SBATCH --time=0-24:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8          # CPUs (threads) per command
#SBATCH --mem=8G                  # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file


module load GCC/11.2.0 GCC/11.3.0 Bowtie/1.3.1

input_dir=/scratch/user/vignesh.bioenv/Metatranscriptomics/FASTQs/All-Samples/Controls
beef_index=/scratch/user/vignesh.bioenv/Metatranscriptomics/Host_genomes/Btaurus

for R1 in $input_dir/*_1.fq.gz
do
#Get the corresponding R2 file
R2=${R1/_1/_2}
bowtie -x $beef_index/genome_index -1 $R1 -2 $R2 --sam $R1.sam
done

