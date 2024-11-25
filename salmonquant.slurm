#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Salmon-Quantification          # job name
#SBATCH --time=7-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=log/stdout.%j          # save stdout to file
#SBATCH --error=log/stderr.%j           # save stderr to file

module load GCC/11.3.0 Salmon/1.10.1

for R1 in *.conrem_1.fq.gz
do
echo "Processing sample ${R1/.conrem_1.fq.gz//}"
#Get the corresponding R2 file
R2=${R1/_1/_2}
salmon quant -i megahit_out_beef/prodigal_out/salmondb/salmon_index -1 $R1 -2 $R2 -p 48 --validateMappings -o megahit_out_beef/prodigal_out/salmon_quant/${R1/.conrem_1.fq.gz//}
done
