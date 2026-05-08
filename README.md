# 🧬 Epigenetic regulation of liver endothelial cells (LSECs) as a novel target to boost immunotherapy efficacy in hepatocellular cancer (HCC)

This project focuses on ***EHMT2***, an epigenetic regulator, and its links with immune pathways and the *endothelial-to-mesenchymal transition* in *liver sinusoidal epithelial cells* (LSECs) which aid in the development of *hepatocellular cancer* (HCC).

## 📋 Actual Work Timeline

```mermaid
gantt
    title Module 6 Project Workflow & Analytical Timeline
    dateFormat  DD-MM-YYYY
    axisFormat  %m/%Y
    tickInterval 1month
    weekday monday

    section Setup & Planning
    Workflow Setup & Comm Channels     :setup, 06-10-2025, 24-10-2025
    Pipeline Selection (nf-core)       :git, 24-10-2025, 03-11-2025

    section Data Analysis
    Raw Data Access & Initial Allocation:milestone, data, 03-11-2025, 1d
    ATAC-Seq Processing                :atac, 03-11-2025, 28-11-2025
    RNA-Seq Task Reassignment          :milestone, reassign, 17-11-2025, 1d
    RNA-Seq Processing                 :rna, 17-11-2025, 12-12-2025
    Multi-Omics Integ.                 :integ, 12-12-2025, 13-01-2026
    HCC scRNA                          :integ, 12-03-2026, 05-04-2026

    section Deliverables
    Video Pitch                        :milestone, video, 15-01-2026, 1d
    Report Drafting 1st Version        :milestone, draft, 28-04-2026, 1d
    Final Report Submission            :milestone, 07-05-2026, 1d
```
## 🎯 Objectives

```mermaid
---
config:
  theme: dark
---
flowchart TB
    A["Raw ATACseq Data"] -- "QC, alignment, peak calling" --> B["Processed ATACseq Data"]
    B-- "Annotation" -->C["Differentially Accessible Regions DARs"]
    I["Raw ATACseq Data"] -- "QC, alignment, quantification" --> D["Processed Bulk RNAseq Data"]
    D["Processed Bulk RNAseq Data"]-- "DESeq2" -->E["Differentially Expressed Genes DEGs"]

    X[ ]:::empty
    C --- X
    E --- X
    X --> F["Candidate Genes\n(DARs linked to DEGs)"]

    Y[ ]:::empty
    F --- Y
    Y --> H["HCC scRNA\npublic dataset"]
    Y --> J["Drug prediction\nMolecular docking"]

    Z[ ]:::empty
    H --- Z
    J --- Z
    Z --> G["Role of EHMT2 in HCC"]

    classDef empty height: 0, width: 0
```

1) Analyse RNA-seq data from LSECs samples with and without *EHMT2* knockdown to identify differentially expressed genes (DEGs).

2) Analyse ATAC-seq data from LSECs samples with and without *EHMT2* knockdown to identify differentially accessible regions (DERs).

3) Compare DEGs against DERs to validate downstream effects of *EHMT2* knockdown.

4) Study DEGs and DERs to find associations.

## 📁 Files & Directories

```bash
bioinformatics-project/
├── ATAC-Seq/                     # ATAC-Seq scriptsand results
├── RNA-Seq/                      # RNA-Seq scriptsand results
├── Integrated_Multiomics/        # Combining ATAC-Seq and RNA-Seq
├── HCC_scRNA/                    # HCC scRNA analysis
├── EnviromentBuilding/           # Including references building
├── Documents/                    # Report draft, slides and meeting records
├── README.md
└── .gitignore                    # list of files to not track
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

