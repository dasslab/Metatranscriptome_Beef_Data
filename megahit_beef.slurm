#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=MEGAHIT-BEEF          # job name
#SBATCH --time=7-00:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360G                  # total memory per node
#SBATCH --output=log/stdout.%j          # save stdout to file
#SBATCH --error=log/stderr.%j           # save stderr to file

module load GCCcore/11.2.0 MEGAHIT/1.2.9

megahit -t 48 -1 B1_8B.conrem_1.fq.gz,T17.conrem_1.fq.gz,T18.conrem_1.fq.gz,T19.conrem_1.fq.gz,T20.conrem_1.fq.gz,T21.conrem_1.fq.gz,T22.conrem_1.fq.gz,T23.conrem_1.fq.gz,T24.conrem_1.fq.gz,T25.conrem_1.fq.gz,T26.conrem_1.fq.gz,T27.conrem_1.fq.gz,T28.conrem_1.fq.gz,T29.conrem_1.fq.gz,T30.conrem_1.fq.gz,T31.conrem_1.fq.gz,T32.conrem_1.fq.gz,T33.conrem_1.fq.gz,T34.conrem_1.fq.gz,T35.conrem_1.fq.gz,T36.conrem_1.fq.gz,T38.conrem_1.fq.gz,T39.conrem_1.fq.gz,T40.conrem_1.fq.gz,T41.conrem_1.fq.gz,T42.conrem_1.fq.gz,T43.conrem_1.fq.gz,T44.conrem_1.fq.gz,T45.conrem_1.fq.gz -2 B1_8B.conrem_2.fq.gz,T17.conrem_2.fq.gz,T18.conrem_2.fq.gz,T19.conrem_2.fq.gz,T20.conrem_2.fq.gz,T21.conrem_2.fq.gz,T22.conrem_2.fq.gz,T23.conrem_2.fq.gz,T24.conrem_2.fq.gz,T25.conrem_2.fq.gz,T26.conrem_2.fq.gz,T27.conrem_2.fq.gz,T28.conrem_2.fq.gz,T29.conrem_2.fq.gz,T30.conrem_2.fq.gz,T31.conrem_2.fq.gz,T32.conrem_2.fq.gz,T33.conrem_2.fq.gz,T34.conrem_2.fq.gz,T35.conrem_2.fq.gz,T36.conrem_2.fq.gz,T38.conrem_2.fq.gz,T39.conrem_2.fq.gz,T40.conrem_2.fq.gz,T41.conrem_2.fq.gz,T42.conrem_2.fq.gz,T43.conrem_2.fq.gz,T44.conrem_2.fq.gz,T45.conrem_2.fq.gz -o megahit_out_beef




