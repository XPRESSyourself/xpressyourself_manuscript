# To Do for BioRXiv
<!--- Finished
Bio insight section
Export single summary images
  - Fix so it doesn't keep subplots as it goes
-->

# Latest specs
- With masking and cufflinks
  - 257 GB input fastq files
  -  u0690617 7130188          rutter 3-00:00:00   06:19:55 2-07:52:58         16    62.50Gn                 0:0
              7130188.bat+     rutter              06:19:55 2-07:52:58         16    62.50Gn  40771596K      0:0
              7130188.ext+     rutter              06:19:58  00:00.003         16    62.50Gn          0      0:0
  - ~ 1.48 min / GB input
  - ~ 11.86 min / input file 

Correlate ATF4 with methianine s35 translation
Not classical reegualtion -- starts out missing ATG, then when translation efficiency poor, can actually get to the ATG, doesnt get knocked off in microORFs
GCN2 has same microORFs?
Plot scatters with distribution to make change points
ATF4 could be repressing or activating, check that
Would they have found same down genes with what they did

Be gracious to original authors, state what they did well and why we chose the dataset, and then tell how methods have improved and show advantages to re-analyzing data with the most current, complete methods

### Before sending manuscript to co-authors
1. Add correlation sample axis labels (Ingolia vs XPRESSpipe)
2. GTFmod test cases
  - Try best to follow Ensembl guidelines
3. Read distribution sub-module
4. Analysis test cases (DESeq2, complexity, read dist, periodicity, metagene)
5. Overall test cases for pipelines
6. Add cufflinks info to paper
  - Run quick test to make sure it flows
  - Cite https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0734-x#CR43

### Before bioRXiv
1. Incorporate batch control wrapper
2. UMI
3. Output example housekeeper with even coverage from bigwig files
  - Make bigwig output standard
4. Test docker instance, ensure R deps are working within
  - Make it upload automatically for each new tag
5. Build documentation
  - Add instructions for beginners
  - Add IGV instructions
6. Make sure aspects of XPRESStools discussed are tested (PCA, volcano plots, normalizations)
7. Solve OOM errors on HPC
8. Figure out analysis output issues on HPC
9. Get test cases from people
10. TCGA validation
11. Add note about truncator only editing exon files, thus must point htseq to exons
  - Could add editing transcripts and genes later
12. Get tests back from Jeff, Alex, Yeyun
13. Incorporate rRNA prober

### After bioRXiv
1. Add instructions on how to add to the project
2. Walkthrough video
  - Add videos like example in this [README](https://github.com/manubot/manubot)
3. More tests!
4. Jupyter notebook examples
5. AWS walkthrough and cost breakdown


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
  - Do lit. research and hint at bio insight, maybe do some validation in drug-course with qPCR and WB, but not much more than that
- NOTCH3 not annotated to be related to ISR from a quick google search


# Test people:
- Jeff Morgan
- Boris Zinshteyn?
- Alex Bott
- Yeyun Ouyang
