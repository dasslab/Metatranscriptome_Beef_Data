#!/bin/bash
#SBATCH --job-name=fastqc           # job name
#SBATCH --time=01:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=2           # CPUs (threads) per command
#SBATCH --mem=12G                   # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file

module purge
module load FastQC/0.11.9-Java-11

fastqc *.fq.gz -o fastqc_raw_reports
