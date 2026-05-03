#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --time=1440:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load Bowtie2/2.5.4-GCC-12.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/References_Index
mkdir Bowtie2_index

bowtie2-build --threads 16 /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/references/GRCh38.p14.genome.fa.gz /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index
