# Title
<b>XPRESSyourself Bioinformatics Suite</b>: Automating and Democratizing High-Throughput Sequencing

# Software Validation
- Beef up test coverage while preparing manuscript
- Have testers run pipeline and analysis notebook and critique documentation
  - XPRESStools
    - Alex Bott
    - Yeyun Ouyang
  - XPRESSpipe
    - Jeff Morgan (riboseq)
    - Mark Wadsworth? (PE RNAseq)

# Biological Validation
- TCGA
  - Take a sample of raw datasets and compare FPKM output with publicly available FPKM values
  - Display one example, and all together linreg scatter plot
  - Have table breaking down each component sample
- Riboprof
  - Run example notebook using this dataset
  - Validate points made in [manuscript](https://www.biorxiv.org/content/biorxiv/early/2018/09/04/408054.full.pdf)
    -
    -
    -
    -
    -
  - Dig out novel components within the metabolism/mitochondria sphere
    - Get collab with iron community from Jared to validate (don't include in manuscript until story is ready, just include validation points until then)
    - Dennis Winge, Jerry Kaplan

# Progress / To Do
- Riboprof truth set
  - Download SRA :ballot_box_with_check:
  - Troubleshoot low RPF alignment rate
  - Quality Control final dataset

- TCGA truth set
  - Get datasets :ballot_box_with_check:
  - Create pipeline Docker :ballot_box_with_check:
  - Implement pipeline in CGC
  - Run alignment in CGC, align to full genome  
  - Compare alignments post counting

- Features
  - Complexity
  - Validate periodicity more
  - Metagene

- Testing
  - Build more test cases
  - Code coverage



# Abstract
With the advent of high-throughput sequencing platforms, expression profiling is becoming commonplace in medical research. However, for the general user, a computational overhead exists. One scenario is the user has no bioinformatics background, and therefore they must pay a bioinformatician to process their data. Another scenario is where the user has bioinformatics experience, but manually process files one by one or are not aware of the proper settings to be used during processing. The XPRESSyourself suite aims to remove these barriers and create a tools to help standardize and increase throughput of data processing and analysis. The XPRESSyourself suite is currently broken down into two software packages. The first, XPRESSpipe, automates the pre-processing, alignment, quantification, normalization, and quality control of single-end and paired-end RNAseq, as well as ribosome profiling sequence data. The second, XPRESStools, is a Python toolkit for expression data analysis, compatible with private or public microarray and RNAseq datasets. This software suite is designed where features can easily be modified, and additional packages can be included for processing of other data types in the future, such as CHIPseq or genome alignment.


# Introduction

### How has RNAseq changed science?
  - Diagnostics ([citation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015936/))
    - Identify CNV, etc undetectable by WGS
    - Expression profiles predict disease sub-type (ie. BRCA)
    - Identify previously unannotated de novo expression sites
    - Aids in diagnosing rare disease and mitochondrial disorders where WES info is limited
  - Identify expression outliers in rare disease ([citation](https://www.biorxiv.org/content/10.1101/408492v1))
    - Often 1 gene contributes to disease
    - Average 1.2 genes per disease sample where z-score > 5


### Need for standardization
  - While some RNAseq methods are mature, still no standard protocol for wet or dry lab aspects
    - [Ingolia 2017 method](https://www.ncbi.nlm.nih.gov/pubmed/28579404), where states 10 yrs in, no standardized protocol
  - Uniform processing allows for better cross-study comparison ([citation](https://www.nature.com/articles/sdata201861))
  - Many workflows rely on antiquated software

### Need for automation
  - [Source 1](https://www.nature.com/articles/s41540-018-0054-3)
    - "transcriptomics data volume has a doubling rate of 7 months"
    - "progressive automation of laboratory work"
    - "a future where data analysis will be the bottleneck in biological research."
  - [Source 2](http://science.sciencemag.org/content/344/6184/653.3/tab-e-letters)
    - "Unprecedented advances have been made in the speed and throughput of next generation sequencing (NGS) platforms"
    - "high demands on the bioinformatics tools necessary for analysis of the data generated, which has grown exponentially"
    - "Although hundreds of thousands of samples have been sequenced, our ability to find, associate, and implicate genetic variants and candidate disease genes far outstrips our ability to understand them"
    - "Many researchers are comfortable with NGS technology, but encounter difficulties with the bioinformatics portion of their workflow"
    - "once clear bioinformatics procedures are established and optimized this bottleneck can be removed, resulting in smooth and routine data interpretation processes and expedited research discoveries"
  - Newer seq platforms outputting even more data at a time
    - Ex: scRNAseq large datasets, need faster methods

### Need for expanded toolkit and what XPRESSyourself offers
  - No standardized software for the following:
    - Periodicity calculation of most abundant fragment
    - Truncate transcript reference to avoid mapping to regions of transcript where known library creation biases exist
    - rRNA depletion -- ability to identify most-abundant or problematic fragments
      - At time of writing, Illumina removing their yeast RiboGold from market as stand-alone kit
      - Commercial kits often insufficient at removing fragments created in ribosome footprinting
  - Some analytical tools not yet available in Python (at least to my knowledge currently)
    - PCA with confidence intervals
    - Automated RPKM/FPKM normalization
    - More automated (single-line) coded options for differential expression analysis and batch effect normalization
    - Set gene threshold (minimum number of reads per gene to keep gene in analysis)
    - Automated volcano plot with different colors for different groups and (semi) auto-labeled labels for outliers
    - Automated linear regression analysis of one gene against all
  - Microarray tools
    - Probe collapse removing multi-mapping probes and averaging multiple probes for a single gene

### Need for a package accessible to all
  - Designed to be flexible, easy to adapt to other sequence techs or add features
  - Thoroughly documented
    - [XPRESSpipe documentation](https://xpresspipe.readthedocs.io/en/latest)
    - [XPRESStools documentation](https://xpresstools.readthedocs.io/en/latest/)
  - Interactive notebooks
    - Jupyter notebook with example analysis (use output as biological insight for paper) that users can tweak to easily run their own analysis

### Competition
  - [ARCHS4](https://www.nature.com/articles/s41467-018-03751-6) and [BioJupies](https://www.ncbi.nlm.nih.gov/pubmed/30447998)
    - Use Kallisto -- less memory but not as accurate
    - Non-flexible
    - Analysis notebooks don't have many customizable arguments
    -
  - Galaxy
    - If using their clusters, data limits (can be navigated by running on own cluster)
    - Many published workflows aren't updated, rely on old software
    - <i><b><u>What are the shortcomings that this package is overcoming?</u></b></i>
    -
    -
  - [Recount](https://f1000research.com/articles/6-1558/v1)
    - Used splice-aware Rail-RNA alignment software
    -
  - [TCGA pipeline tool](https://github.com/akahles/icgc_rnaseq_align)
    - Limited functionality
    - No analysis
    - Personally couldn't get it to run out of the box
    -


# Outline

### [XPRESSpipe repository](https://github.com/XPRESSyourself/XPRESSpipe) and [documentation](https://xpresspipe.readthedocs.io/en/latest)
  - <b>Figure 1</b> overviewing inputs and outputs of pipeline
  ![XPRESSpipe](https://raw.githubusercontent.com/XPRESSyourself/XPRESSpipe/master/docs/content/xpresspipe_overview.png)
  - Run each step of processing independently
    - Trim
      - FastP
    - Align
      - 2-pass STAR
    - Count
      - HTSeq-Count union exon method
    - Normalize/Scale
      - Batch
      - RPM, RPKM, FPKM
    - Metagene/Quality Control
      - Read distributions
      - Metagene profiles
      - 3-nt periodicity
      - MultiQC
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

### [XPRESStools repository](https://github.com/XPRESSyourself/XPRESStools) and [documentation](https://xpresstools.readthedocs.io/en/latest/)
  - Add workflow figures, see example outputs in documentation
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
  - Parallelize dataframe computation

### Speed and Cost
  - <b>Figure 3</b>: Cost and speed breakdown?
  - Calculate cost on Amazon cluster
  - Calculate speed per 1 million reads for the pipeline


# Potential Biological Insights Section
  - <b>Figure 4</b>
  - Run several GEO datasets through pipeline and analysis and pull out some novel insight off topic from original papers
    - Find datasets that could complement each other and parallel process


# Methods

### [XPRESSpipe](https://github.com/XPRESSyourself/XPRESSpipe)
  - Normalization functions
  - Multiprocessing optimization
  - Dependencies and logic behind using each package
  - Default parameters for STAR, etc
  - Calculating most abundant read length for periodicity
  - rRNA probe design
  -
  -
  -

### [XPRESStools](https://github.com/XPRESSyourself/XPRESStools)
  - Confidence intervals
  - Collapser
  - Normalization functions
  - Automated 1vsAll linear regression modeling
  -
  -
  -






# Notes

### INTRO
RNAseq and riboseq powerful tools   

Automated analysis critical   
Pipelines that exist that tackle this   

Few options for ribosome profiling   
Talk about additional tools    

Describe a new protocol for both, plus new tools   

Can focus on ribosome profiling and mention also works with RNAseq   
Brought it all together in an automated way   
Flexible easy to use   

Talk about new stuff first in results   

Download someone elses data and find new things, find other references to help explain if possible   
Find same things at a minimum   
Find a metabolism slant   

Don't make RNAseq the focus   

Output stdout and error into log file    
Timestamp or experiment name   
Scan each output for errors and stop if true   

### Quality control:
Measure of complexity   
	- Report from alignment -- # of reads to # of unique locations
	- Saturation calculation
	- When big differences between libraries
	- Random sample a million reads and see how many uniquely align

Multiplexing   
UMI -- later, may be more effort, no standard protocol   

One or two confirming examples for biorxiv   
Codecov and docs   
Get a couple people to test   

### Benchmarking -- not needed for biorxiv
With new capabilities, not really clear who to benchmark against   
	- Focus on new tools
	- Take a couple TCGA samples and run linreg
