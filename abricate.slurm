#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=ABRICATE-RUN-PRODIGAL          # job name
#SBATCH --time=1-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=log/stdout.%j          # save stdout to file
#SBATCH --error=log/stderr.%j           # save stderr to file

module load GCC/8.3.0 OpenMPI/3.1.4 ABRicate/0.9.9

abricate proteins.l500.contigs.fasta --csv abricate.out.csv --threads 8 