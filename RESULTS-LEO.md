# ATAC-Seq and RNA-Seq Integrated Analysis Workflow: Full Document Explanation

This document provides a comprehensive explanation of the multi-omics
workflow illustrated in the original slides (Process.pdf). It describes
how ATAC-Seq and RNA-Seq datasets were processed and integrated to
identify putative regulatory genes associated with EHMT2 knockdown (KD).
The explanation is organized page-by-page according to the source
document.

------------------------------------------------------------------------

## **Page 1 --- Overview of Sample Grouping and ATAC-Seq Processing**

### **Sample Grouping and Key Comparisons**

The top-left section summarizes the experimental design, including
control and KD groups, each with multiple biological replicates (B and
C).\
The red-highlighted comparisons indicate the specific contrasts used for
downstream ATAC-Seq differential accessibility analysis. These
comparisons form the basis for identifying KD-induced chromatin
accessibility changes.

### **ATAC-Seq Analytical Workflow**

The lower portion of Page 1 outlines the pipeline for identifying
putative target genes using two analytical strategies.

### **Method 1: Annotation First, Then Intersection**

1.  MACS2 is used for peak calling on each biological replicate.
2.  Peaks from each replicate are directly annotated to the human genome
    (GRCh38.p14).
3.  Differentially Accessible Areas (DAAs) are defined by intersecting
    annotated peaks between replicates B and C.
4.  Genes linked to these intersected peaks are identified as
    DAA-associated genes.

This method typically yields a larger number of genes.

### **Method 2: Intersection First, Then Annotation**

1.  MACS2 peak calling is performed on each replicate.
2.  Shared peak regions between replicates B and C are identified.
3.  These consensus peak regions are then annotated.
4.  Genes associated with these regions are extracted as DAA-associated
    genes.

This method identifies fewer but more reproducible peaks across
replicates.

In summary, **Method 1 is more inclusive**, while **Method 2 is more
stringent**.

------------------------------------------------------------------------

## **Page 2 --- DAA Identification: Open vs. Closed Chromatin Regions**

Two Venn diagrams compare Method 1 and Method 2 results.

### **Figure A --- Regions Open After KD (KD vs. Control)**

These regions exhibit **increased chromatin accessibility** following
EHMT2 KD.\
The overlap between methods indicates that many regions are consistently
detected, though Method 1 has greater sensitivity.

### **Figure B --- Regions Closed After KD (Control vs. KD)**

These represent genomic regions with **reduced accessibility** after
EHMT2 KD.\
Again, Method 1 produces a larger gene set, whereas Method 2 yields a
refined subset of reproducible accessibility losses.

------------------------------------------------------------------------

## **Page 3 --- Multi-Omics Integration: Open Chromatin and Upregulated Genes (OpenUP)**

Page 3 integrates ATAC-Seq open-region results with RNA-Seq differential
expression from KD1 and KD2.

*RNA-Seq significance thresholds:*\
- \|log₂FC\| \> 1\
- p \< 0.05

### **Figure A --- Intersection Analysis**

-   **ATAC-Seq:** Intersection of open-region results from Methods 1 and
    2\
-   **RNA-Seq:** Intersection of significantly upregulated genes in KD1
    and KD2

This identifies genes that show **both structural chromatin opening and
transcriptional upregulation**, representing strong candidates activated
by KD.

### **Figure B --- Union Analysis**

-   **ATAC-Seq:** Union of open regions from both methods\
-   **RNA-Seq:** Union of upregulated DEGs from both KD experiments

This produces a broader set of candidate genes, capturing genes
supported by either dataset.

------------------------------------------------------------------------

## **Page 4 --- Multi-Omics Integration: Closed Chromatin and Downregulated Genes (CloseDOWN)**

This page mirrors the logic of Page 3 but focuses on genes whose
chromatin and expression are both repressed.

### **Figure A --- Intersection Analysis**

-   **ATAC-Seq:** Intersection of Method 1 and Method 2 for closed
    regions\
-   **RNA-Seq:** Intersection of significantly downregulated genes in
    KD1 and KD2

These genes show consistent reductions in both chromatin accessibility
and expression.

### **Figure B --- Union Analysis**

-   **ATAC-Seq:** Union of closed regions from both methods\
-   **RNA-Seq:** Union of downregulated DEGs

This yields a more inclusive list of potentially EHMT2-regulated
repressed genes.

------------------------------------------------------------------------

## **Page 5 --- Pearson Correlation Heatmaps with EHMT2**

The final page displays two Pearson correlation heatmaps constructed
using RNA-Seq expression values from an EHMT2 tsv file.

### **Figure A --- Genes Open in ATAC-Seq and Upregulated in RNA-Seq**

This heatmap shows correlations between **EHMT2** and the **11 genes**
identified in the union analysis on Page 3B.\
The correlation patterns help infer potential regulatory
interactions---positive or negative---associated with EHMT2 KD.

### **Figure B --- Genes Closed in ATAC-Seq and Downregulated in RNA-Seq**

This heatmap presents correlations between EHMT2 and the **8 genes**
identified in Page 4B.\
These genes may reflect EHMT2-dependent chromatin silencing effects.

------------------------------------------------------------------------

# **Overall Summary**

Across all pages, the workflow integrates structural (ATAC-Seq) and
transcriptional (RNA-Seq) changes to systematically identify
EHMT2-related regulatory genes. Key features include:

1.  Two complementary ATAC-Seq pipelines providing sensitive versus
    stringent detection of differential chromatin accessibility.\
2.  Integration with RNA-Seq to identify genes coherently regulated at
    both chromatin and expression levels.\
3.  Pearson correlation analysis to evaluate co-regulation with EHMT2
    and strengthen candidate gene prioritization.

This multi-omics strategy enables robust identification of potential
EHMT2-regulated genes and biological pathways influenced by knockdown
effects.
