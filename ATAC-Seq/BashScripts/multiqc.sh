#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load MultiQC/1.9-foss-2019b-Python-3.7.4

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC

multiqc ./fastqc

