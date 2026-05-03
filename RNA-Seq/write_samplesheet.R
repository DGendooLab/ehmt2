library(tidyverse)

path = "/rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/fastq/"

samples <- c("A", "B", "C")

conditions <- c(1:6)

samplesheet <- data.frame(
  sample = character(),
  fastq_1 = character(),
  fastq_2 = character(),
  strandedness = character()
)

for (sample_id in samples){
  for (condition in conditions){
      sample_condition <- paste0(sample_id,"_",condition)
      R1 <- paste0(path, "SR-",  sample_id, condition, "_R1_001.fastq.gz")
      R2 <- paste0(path, "SR-",  sample_id, condition, "_R2_001.fastq.gz")
      
      samplesheet <- samplesheet %>%
        add_row(sample = sample_condition,
                fastq_1 = R1,
                fastq_2 = R2,
                strandedness = "auto"
                )
  }
}

write.csv(samplesheet, file = "./samplesheet.csv", row.names = F)
quit()
