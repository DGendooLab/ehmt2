#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --time=1440:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2022a
module load deepTools/3.5.2-foss-2022a
module load SAMtools/1.16.1-GCC-11.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/

# samtools index -@ 16 B1-1-C_rmdup.bam
# samtools index -@ 16 B2-1-KD_rmdup.bam
# samtools index -@ 16 B3-1-C-T_rmdup.bam
# samtools index -@ 16 B4-1-KD-T_rmdup.bam
# samtools index -@ 16 C1-1-C_rmdup.bam
# samtools index -@ 16 C2-1-KD_rmdup.bam
# samtools index -@ 16 C3-1-C-T_rmdup.bam
# samtools index -@ 16 C4-1-KD-T_rmdup.bam

bamCoverage --bam B1-1-C_rmdup.bam -o B1-1-C_rmdup.bw --binSize 10
bamCoverage --bam B2-1-KD_rmdup.bam -o B2-1-KD_rmdup.bw --binSize 10
bamCoverage --bam B3-1-C-T_rmdup.bam -o B3-1-C-T_rmdup.bw --binSize 10
bamCoverage --bam B4-1-KD-T_rmdup.bam -o B4-1-KD-T_rmdup.bw --binSize 10
bamCoverage --bam C1-1-C_rmdup.bam -o C1-1-C_rmdup.bw --binSize 10
bamCoverage --bam C2-1-KD_rmdup.bam -o C2-1-KD_rmdup.bw --binSize 10
bamCoverage --bam C3-1-C-T_rmdup.bam -o C3-1-C-T_rmdup.bw --binSize 10
bamCoverage --bam C4-1-KD-T_rmdup.bam -o C4-1-KD-T_rmdup.bw --binSize 10