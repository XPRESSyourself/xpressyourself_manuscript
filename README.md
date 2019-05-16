# To Do for BioRXiv
<!--- Finished
Bio insight section
Export single summary images
  - Fix so it doesn't keep subplots as it goes
-->

### Before bioRXiv
1. Add correlation sample axis labels (Ingolia vs XPRESSpipe)
2. GTFmod test cases
  - Try best to follow Ensembl guidelines
3. Read distribution sub-module
4. Analysis test cases (DESeq2, complexity, read dist, periodicity, metagene)
5. Overall test cases for pipelines
6. UMI
7. Incorporate batch control wrapper
9. Output example housekeeper with even coverage from bigwig files
  - Make bigwig output standard
10. Test docker instance, ensure R deps are working within
  - Make it upload automatically for each new tag
11. Build documentation
  - Add instructions for beginners
  - Add IGV instructions
12. Make sure aspects of XPRESStools discussed are tested (PCA, volcano plots, normalizations)
13. Solve OOM errors on HPC
14. Figure out analysis output issues on HPC


### After bioRXiv
1. Add instructions on how to add to the project
2. Walkthrough video
  - Add videos like example in this [README](https://github.com/manubot/manubot)
3. Biology validation
4. TCGA validation
5. More tests!
6. Jupyter notebook examples
7. AWS walkthrough and cost breakdown




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
