# To Do for BioRXiv
- Make own read distribution sub-module
- UMI module
- Add FastQC to rRNAprobe
- Find batch control best package
- Export single summary images
  - Fix so it doesn't keep subplots as it goes
- For longest transcript, measure in exon space
- (COMPLETE) Fix deduplication and use deduplicated file for counting, meta-analysis, periodicity
- (COMPLETE) For meta position penalty, only count if in exon space
- Test cases for pipeline
- Test cases for analysis
- Test docker instance, ensure R deps are working within
- Build documentation
  - Add instructions for beginners
  - Add instructions on how to add to the project
- Dataset test case
- Get reference formatted fixed and make in-text superscript


# Citation Notes
- DOI:10.1038/s41598-017-01617-3
  - Benchmarking of RNA-sequencing analysis workflows using whole transcriptome RT-qPCR expression data
  - All quantification methods about similar benchmarked to qPCR validation

# Other Notes
- Single-pass is about 15 min per file
- Two-pass is about 55 min per file
- Plans:
  - Focus on software and new tools
  - When doing bio insight, be clear that this is a tool for showcasing pipeline utility
  - Do lit research and hint at bio insight, maybe do some validation in drug-course with qPCR and WB, but not much more than that

# Title
<b>XPRESSyourself Bioinformatics Suite</b>: Automating and Democratizing High-Throughput Sequencing

# Intro


# Methods
### Reference Curation
- Methods, flexibility, assumptions
- How GTF is modified

### Trim
- Methods, flexibility, assumptions

### Align
- Methods, flexibility, assumptions

### Count
- Methods, flexibility, assumptions
- Swap for RSEM

### Quality Control
##### Read distribution summary

##### Complexity summary
- Uses dupRadar on each file and plot together

##### Metagene summary
- For each read, take left-most coordinate + 1/2 * mapped read length
- Find distance from strand-aware 5' end of transcript range it falls into and take percentage position
- Penalty for multimappers
- Custom refFlat

##### Periodicity summary
- For each read, take position of 3' end mapped - 16 (refs)
- Find distance from strand-aware 5' end of transcript range it falls if first 200 nt from start
- Penalty for multimappers
- Custom refFlat

### Other features
##### rRNAprobe
- Take FastQC output zip files and parse for footprint samples and find over-represented sequences
  - How to handle with no name saying that

##### Normalization and Batch Effect
- RPM, TPM, RPKM, FPKM, Log
  - Automated with GTF used in mapping
- Batch effect

##### DESeq2
- Python wrapper for DESeq2

### Further analysis with XPRESStools
##### Unique Features
- Python volcano plot, PCA confidence plotting
- Ease of plotting uses

##### How to use
- Jupyter notebook example
- Show how to take dataframe dictionary to remap sample name with meaningful name


# Validation
### Replicating publicly available datasets
- Replicate TCGA FPKM counts
- Replicate Ingolia ISRIB plots

### Case Study: Ease of Use
- Questionaire
  - Explain background
  - Explain time catching up with codecademy
  - Explain time to get using
  - Compare output vs actual, were they able to extract key points

### Cost
- AWS example

### Documentation
- Walkthroughs for beginners
  - Add videos like example in this [README](https://github.com/manubot/manubot)
- Explanation on how to build upon
- Test examples for all functions

# Biological Insight
- Do after BioRXiv submission


# Conclusions


# References
- Cite all software used to create this (all dependencies)



# Test people:
- Jeff Morgan
- Boris Zinshteyn?
- Alex Bott
- Yeyun Ouyang
- Jon Belyeu
- Mike Howard
- Jay Gertz
- Cameron Waller
- Others in Dept?
