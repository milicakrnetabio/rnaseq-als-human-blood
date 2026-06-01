# rnaseq-als-human-blood
RNA-seq Analysis: Sporadic ALS vs. Healthy Controls
Whole Peripheral Blood Transcriptomics — GSE234297
� � � � �
Overview
Differential gene expression analysis of whole peripheral blood from sporadic ALS (sALS) patients compared to healthy controls. This is a translational study — unlike mouse spinal cord models, whole blood is a clinically accessible, non-invasive sample type directly applicable to human biomarker discovery.
Pipeline: R (DESeq2) + Python (matplotlib/seaborn)
Biological Background
Sporadic ALS accounts for ~90% of all ALS cases, with no identified single genetic cause. Whole blood transcriptomics offers a window into systemic disease processes — immune dysregulation, inflammatory signaling, and RNA metabolism — that may reflect or drive motor neuron degeneration.
Key pathways identified in this dataset (Maksour et al., 2024):
Ferroptosis activation
Immune and inflammatory dysregulation
Altered leukocyte proportions
Differential transcript usage
Key question: Can whole blood gene expression distinguish sALS patients from healthy controls, and which biological pathways are most affected?
Dataset
Field
Details
GEO Accession
GSE234297
Organism
Homo sapiens
Tissue
Whole peripheral blood (PAXgene RNA tubes)
Groups
sALS (n=96) vs. Healthy controls (n=48)
Sequencing
Illumina short-read, paired-end, ~50M reads/sample
Library prep
Stranded total RNA (RiboZero Plus)
Reference genome
GRCh38 / hg38
Published
Annals of Clinical and Translational Neurology, 2024
How This Differs from the Mouse Project
Feature
Mouse project (GSE38820)
This project (GSE234297)
Organism
Mus musculus
Homo sapiens
Sample type
Spinal cord (LCM)
Whole blood
Invasiveness
Terminal, invasive
Non-invasive
Disease stage
Presymptomatic
Symptomatic
Sample size
n=2/group
n=96 vs. 48
Study type
Pretklinički
Translacijski
This project represents the translational step — asking whether molecular signatures found in mouse models are also present in human patients, using a clinically accessible sample.
Project Structure
rnaseq-als-human-blood/
│
├── R/
│   └── 01_deseq2_analysis.R      # DESeq2 pipeline — runs first
│
├── Python/
│   └── 02_visualization.py       # Visualization — runs second
│
├── data/
│   └── GSE234297_counts.txt      # Download from GEO (see Usage)
│
├── results/                      # Generated on run
│   ├── deseq2_results.csv
│   ├── significant_DEGs.csv
│   ├── normalized_counts_vst.csv
│   ├── sample_metadata.csv
│   ├── plot_pca.png
│   ├── plot_volcano.png
│   ├── plot_heatmap.png
│   ├── plot_summary.png
│   └── session_info.txt
│
└── README.md
Pipeline
┌─────────────────────────────────┐
                    │   R / 01_deseq2_analysis.R       │
                    ├─────────────────────────────────┤
                    │ 1. GEO download / local import   │
                    │ 2. Pre-filter (≥10 counts, ≥48s) │
                    │ 3. DESeq2                         │
                    │    ├─ Size factors                │
                    │    ├─ Dispersion estimation       │
                    │    ├─ Negative binomial GLM       │
                    │    └─ Wald test + BH-FDR          │
                    │ 4. LFC shrinkage (apeglm)         │
                    │ 5. VST normalization              │
                    │ 6. Export CSVs → results/         │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │   Python / 02_visualization.py   │
                    ├─────────────────────────────────┤
                    │ 1. PCA (top 500 variable genes)  │
                    │ 2. Volcano plot                   │
                    │ 3. Heatmap (top 50 DEGs)          │
                    │ 4. Summary panel                  │
                    └─────────────────────────────────┘
Usage
1. Clone repository
git clone https://github.com/YOUR_USERNAME/rnaseq-als-human-blood.git
cd rnaseq-als-human-blood
2. Download count matrix
1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234297
2. Scroll to "Supplementary file"
3. Download count matrix
4. Save as: data/GSE234297_counts.txt
Or let R download automatically:
source("R/01_deseq2_analysis.R")  # will attempt GEOquery download
3. R dependencies
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery", "apeglm"))
install.packages(c("dplyr", "tibble", "readr", "stringr"))
4. Run R analysis
setwd("rnaseq-als-human-blood/")
source("R/01_deseq2_analysis.R")
5. Python dependencies
pip install pandas numpy matplotlib seaborn scikit-learn
6. Run Python visualization
python Python/02_visualization.py
Statistical Notes
Parameter
Value
Normalization
Median-of-ratios (DESeq2)
Test
Wald test (negative binomial GLM)
LFC shrinkage
apeglm
FDR correction
Benjamini-Hochberg
Significance
padj < 0.05, |log2FC| > 0.58 (≈1.5x)
Pre-filter
≥10 counts in ≥48 samples
VST
For PCA and heatmap only
Why |LFC| > 0.58? 0.58 ≈ log2(1.5) — corresponds to 1.5-fold change. With large n (n=144 total), DESeq2 has high power to detect small but biologically irrelevant changes. The LFC threshold filters for biologically meaningful differences.
Candidate Biological Pathways
Based on ALS blood transcriptomics literature:
Category
Genes
Ferroptosis
GPX4, SLC7A11, ACSL4, PTGS2, HMOX1
Immune / Inflammatory
S100A8, S100A9, S100A12, LYZ, IL6R
Complement
C1QB, C3, CFB
RNA metabolism
TARDBP, FUS, HNRNPA1, STMN2
Neurofilaments
NEFL, NEFM, NEFH
Mitochondrial
TFAM, PINK1, PARK7
Reproducibility
Resource
Location
Raw data
NCBI GEO GSE234297
Code
This repository (GitHub)
R versions
results/session_info.txt
Parameters
This README
Ontologies
UBERON:0000178 (blood), DOID:332 (ALS), NCBITaxon:9606
This project follows FAIR principles — data is Findable (GEO accession + DOI), Accessible (public), Interoperable (standard formats, ontologies), and Reusable (documented parameters, open license).
References
Maksour S et al. (2024). RNA sequencing of peripheral blood in amyotrophic lateral sclerosis reveals distinct molecular subtypes. Annals of Clinical and Translational Neurology. PMC10946588
Love MI, Huber W, Anders S (2014). DESeq2. Genome Biology. doi:10.1186/s13059-014-0550-8
Zhu A et al. (2019). apeglm. Bioinformatics. doi:10.1093/bioinformatics/bty895
Author
Milica Krneta — Biologist / Data Curator MSc Mycology & Microbiology, University of Belgrade Skills: RNA-seq, DESeq2, data curation, FAIR principles, Python, R
Part of a genomics data curation portfolio. See also:
ALS Motor Neurons — Mouse Model (GSE38820)
COVID-19 vs. Healthy Controls (GSE157103)