### LaTeX manuscript files for XPRESSyourself


### Abstract
With the advent of high-throughput sequencing platforms, expression profiling is becoming commonplace in medical research. However, for the general user, a computational overhead exists. One scenario is the user has no bioinformatics background, and therefore they must pay a bioinformatician to process their data. Another scenario is where the user has bioinformatics experience, but manually process files one by one or are not aware of the proper settings to be used during processing. The XPRESSyourself suite aims to remove these barriers and create a tools to help standardize and increase throughput of data processing and analysis. The XPRESSyourself suite is currently broken down into two software packages. The first, XPRESSpipe, automates the pre-processing, alignment, quantification, normalization, and quality control of single-end and paired-end RNAseq, as well as ribosome profiling sequence data. The second, XPRESStools, is a Python toolkit for expression data analysis, compatible with private or public microarray and RNAseq datasets. This software suite is designed where features can easily be modified, and additional packages can be included for processing of other data types in the future, such as CHIPseq or genome alignment.


### Introduction

# Need for standardization
  - Reference [Ingolia 2017 method](https://www.ncbi.nlm.nih.gov/pubmed/28579404), where states 10 yrs in, no standardized protocol
  - Uniform processing allows for better cross-study comparison

# Need for automation
  - ... there is an exponentially increasing amount of functional genomics data… transcriptomics data volume has a doubling rate of 7 months, and high-throughput workflows for proteomics and metabolomics are becoming increasingly available. Furthermore, the miniaturization of these techniques and the progressive automation of laboratory work through microfluidics chips promises a future where data analysis will be the bottleneck in biological research. (Costello & Martin, NPJ Systems Biology and Applications, 2018)
  - Unprecedented advances have been made in the speed and throughput of next generation sequencing (NGS) platforms over the last decade. This progress has imposed increasingly high demands on the bioinformatics tools necessary for analysis of the data generated, which has grown exponentially. Although hundreds of thousands of samples have been sequenced, our ability to find, associate, and implicate genetic variants and candidate disease genes far outstrips our ability to understand them. Many researchers are comfortable with NGS technology, but encounter difficulties with the bioinformatics portion of their workflow, rendering NGS a less attractive option as their primary sequencing platform. However, once clear bioinformatics procedures are established and optimized this bottleneck can be removed, resulting in smooth and routine data interpretation processes and expedited research discoveries. (Funari and Canosa, Science Webinar Summary, 2014)

# Need for additional tools (ie. ribosome profiling)
  - No standardized software for the following:
    - Periodicity calculation of most abundant fragment
    - Truncate transcript reference to avoid mapping to regions of transcript where known library creation biases exist
    - rRNA depletion -- ability to identify most-abundant or problematic fragments
      - At time of writing, Illumina removing their yeast RiboGold from market as stand-alone kit
      - Commercial kits often insufficient at removing fragments created in ribosome footprinting

# Need for a package accessible to all
  - Designed to be flexible, easy to adapt to other sequence techs
  - Easy to add
  - Thoroughly documented

# Competition
  - [ARCHS4](https://www.nature.com/articles/s41467-018-03751-6) and [BioJupies](https://www.ncbi.nlm.nih.gov/pubmed/30447998)
    - Kallisto
    - non-flexible?  
  - Galaxy -- its Galaxy
  - [Recount] (https://f1000research.com/articles/6-1558/v1)
  - [TCGA pipeline tool](https://github.com/akahles/icgc_rnaseq_align)
    - Limited functionality
    - No analysis
    - Buggy


### Outline

# [XPRESSpipe](https://github.com/XPRESSyourself/XPRESSpipe)
  - <b>Figure 1</b>: Pipeline and output overview
  - Run each step of processing independently
    - Trim
    - Align
    - Count
    - Normalize
    - Metagene/Quality Control
  - Run as a pipeline in an automated, parallel manner
    - Single-end RNAseq
    - Paired-end RNAseq
    - Ribosome Profiling
  - STAR
    - [TCGA standard](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
    - [Best output quality](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/)
    - Map splice junctions for better alignment
  - Automates
    - Normalization
    - Differential expression
    - Meta-analysis
      - Read size distributions
      - Periodicity of dominant footprint length
      - Meta-gene analysis
      - Pipeline will summarize steps using MultiQC
    - Reference curation
      - STAR reference
      - refFLAT
      - Truncated or coding-only quantification reference files
    - Identification of over-abundant rRNA fragments for depletion

# [XPRESStools](https://github.com/XPRESSyourself/XPRESStools)
  - <b>Figure 2</b>: Test data and example outputs
  - Access data
  - Format for processing
  - Normalize for downstream analysis
    - Probe collapse for microarray datasets
  - Perform quality control
  - Analyze
    - Clustered heatmaps of all genes or subset
    - Single- or multi-gene analysis for sample types
    - Scatterplots and linear regressions
    - Volcano plot
      - Label genes, plot thresholds, highlight gene sets
    - PCA
      - 2D or 3D
      - Plot confidence intervals
      - Interactive optional

# Speed and Cost
  - <b>Figure 3</b>: Cost and speed breakdown
  - Calculate cost on Amazon cluster
  - Calculate speed per 1 million reads for the pipeline


### Methods

# XPRESSpipe
  - Normalization functions
  - Multiprocessing optimization
  -

# XPRESStools
  - Confidence intervals
  - Collapser
  - Normalization
  -
