# 🧬 Epigenetic regulation of liver endothelial cells (LSECs) as a novel target to boost immunotherapy efficacy in hepatocellular cancer (HCC)

This project focuses on **EHMT2**, an epigenetic regulator, and its links with immune pathways and the *endothelial-to-mesenchymal transition* in *liver sinusoidal epithelial cells* (LSECs) which aid in the development of *hepatocellular cancer* (HCC).

## 📋 Plan Overview

```mermaid
gantt
    title Semester 1 Timeline
    dateFormat  DD-MM-YYYY
    axisFormat  %d/%m
    tickInterval 1week
    weekday monday

    section Video & Report
    Literature Review + Workflow Setup       :lit, 09-10-2025, 22-10-2025
    1st Meeting w/ Deena     :milestone, deena, 22-10-2025, 1d
    Video Pitch              :pitch, after real, 13-12-2025
    Video Pitch Deadline             :milestone, 15-01-2026, 1d
    Report 1st Draft         :after pitch, 12-01-2026

    section Analysis
    Real DAA + DEA           :real, after deena, 5w
    Visualisations    :after real, 13-12-2025

    section Publication Push
    Practise DEA + DAA                 :practise, 13-10-2025, 22-10-2025
    Propose Publication w/ Deena      :proposal, after practise, 3w
    Scout Public DB          :db, after pitch, 1w
    Meta Analysis            :meta, after db, 12-01-2026
```
## 🎯 Objectives

```mermaid
---
config:
  theme: dark
---
flowchart TB
    A["Raw ATACseq Data"] -- "QC, peak calling" --> B["Processed ATACseq Data"]
    B-- "DAA w/ DiffBind" -->C["Differentially Accessible Regions DARs"]
    D["Processed Bulk RNAseq Data"]-- "DEA w/ DESeq2" -->E["Differentially Expressed Genes DEGs"]

    X[ ]:::empty
    C --- X
    E --- X
    X --> F["DARs linked to DEGs"]

    Y[ ]:::empty
    F --- Y
    Y --> G["Role in EndoMT and immune pathways"]
    H["HCC Metadata Immune Markers"] --- Y

    classDef empty height: 0, width: 0
```

1) Analyse RNA-seq data from LSECs samples with and without EHMT2 knockdown to identify differentially expressed genes (DEGs).

2) Analyse ATAC-seq data from LSECs samples with and without EHMT2 knockdown to identify differentially accessible regions (DERs).

3) Compare DEGs against DERs to validate downstream effects of EHMT2 knockdown.

4) Study DEGs and DERs to find associations with immune markers.

## 📁 Files & Directories

```bash
bioinformatics-project/
├── ATAC-Seq/                     # ATAC-Seq scriptsand results
├── RNA-Seq/                      # RNA-Seq scriptsand results
├── Integrated_Multiomics/        # Combining ATAC-Seq and RNA-Seq
├── HCC_scRNA/                    # HCC scRNA analysis
├── EnviromentBuilding/           # Including references building
├── Documents/                    # Report draft, slides and meeting records
└── README.md
└── .gitignore                # list of files to not track
```
## ⚠️ Notice
Cause of limitation of file size in GitHub, files listed below are not included in repo.

- .fastq.gz (Raw Data)
- .RData (Image of analysis)
- .RDS Larger than 100 MB
- .bedGraph (MACS2 visualizaion data of ATAC-Seq)
- .fasta, .gtf, .gff3 (Genome references files)



## 👥 Contributors

- Zhaoshuo Liu
- Miguel Alburo
- Simran Panda
- Yash Dhiman

