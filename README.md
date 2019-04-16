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


# To Do for BioRXiv
- Fix deduplication and use deduplicated file for counting, meta-analysis, periodicity
- For meta position penalty, only count if in exon space
- Test cases for pipeline
- Test cases for analysis
- Test docker instance
- Build documentation
  - Add instructions for beginners
  - Add instructions on how to add to the project
- Dataset test case
- Get reference formatted fixed and make in-text superscript
- Implement UMI to trim
