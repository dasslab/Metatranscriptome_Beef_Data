#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=BEEF_CONTAMINATION_REMOVAL_CONTROLS
#SBATCH --time=3-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=stdout.%x.%j
#SBATCH --error=stderr.%x.%j

module load GCC/11.2.0 GCC/11.3.0 Bowtie/1.3.1

control_index=/scratch/user/vignesh.bioenv/Metatranscriptomics/FASTQs/All-Samples/Controls/megahit_out_controls/control_index

input_dir=/scratch/user/vignesh.bioenv/Metatranscriptomics/FASTQs/All-Samples/Beef/Host_removed_reads/GZIP

for R1 in $input_dir/*_1.unaligned.fq.gz
do
#Get the corresponding R2 file
R2=${R1/_1/_2}
bowtie -x $control_index/genome_index -1 $R1 -2 $R2 --un controls_mapped/$R1.conrem.fq
done





