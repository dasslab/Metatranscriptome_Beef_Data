#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=Kraken2           # job name
#SBATCH --time=3-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=Kraken2out/hprclog/stdout.%j          # save stdout to file
#SBATCH --error=Kraken2out/hprclog/stderr.%j           # save stderr to file


module load GCC/10.3.0 OpenMPI/4.1.1 Kraken2/2.1.2

for R1 in *.fq.gz
do 
R2=${R1/_1/_2}
kraken2 --threads 96 --db Kraken2/standard --output Kraken2out/${R1/_1.fq.gz/}_output.txt --report Kraken2out/report/${R1/.conrem_1.fq.gz/}.report.txt --paired $R1 $R2
done