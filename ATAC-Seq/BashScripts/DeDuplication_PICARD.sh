#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --time=4320:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load SAMtools/1.21-GCC-12.3.0
module load picard/3.4.0-Java-17

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/

# samtools addreplacerg -r "@RG\tID:B1-1-C_sorted\tSM:B1-1-C\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B1-1-C_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B1-1-C_sorted.bam
# java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B1-1-C_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.log --REMOVE_DUPLICATES true 

# samtools addreplacerg -r "@RG\tID:B2-1-KD_sorted\tSM:B2-1-KD\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B2-1-KD_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B2-1-KD_sorted.bam
# java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B2-1-KD_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.log --REMOVE_DUPLICATES true 

# samtools addreplacerg -r "@RG\tID:B3-1-C-T_sorted\tSM:B3-1-C-T\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B3-1-C-T_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B3-1-C-T_sorted.bam
# java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B3-1-C-T_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.log --REMOVE_DUPLICATES true 

# samtools addreplacerg -r "@RG\tID:B4-1-KD-T_sorted\tSM:B4-1-KD-T\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B4-1-KD-T_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B4-1-KD-T_sorted.bam
# java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B4-1-KD-T_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.log --REMOVE_DUPLICATES true 

samtools addreplacerg -r "@RG\tID:C1-1-C_sorted\tSM:C1-1-C\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C1-1-C_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C1-1-C_sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C1-1-C_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.log --REMOVE_DUPLICATES true 

samtools addreplacerg -r "@RG\tID:C2-1-KD_sorted\tSM:C2-1-KD\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C2-1-KD_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C2-1-KD_sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C2-1-KD_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.log --REMOVE_DUPLICATES true 

samtools addreplacerg -r "@RG\tID:C3-1-C-T_sorted\tSM:C3-1-C-T\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C3-1-C-T_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C3-1-C-T_sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C3-1-C-T_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.log --REMOVE_DUPLICATES true 

samtools addreplacerg -r "@RG\tID:C4-1-KD-T_sorted\tSM:C4-1-KD-T\tPL:Illumina\tLB:Hg38" -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C4-1-KD-T_sortedRg.bam /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C4-1-KD-T_sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C4-1-KD-T_sortedRg.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.bam -M /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.log --REMOVE_DUPLICATES true 