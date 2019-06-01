# To Do for BioRXiv
<!--- Finished
Bio insight section
Export single summary images
  - Fix so it doesn't keep subplots as it goes
-->

### Before bioRXiv
Week ending June 2:
[X] Genome size in curateReference
[X] Set Cufflinks default
[?] Dedup with fixmate step before, not default
  - Apparently not needed with coordinate sorted
[ ] Longest transcript over exon space
[ ] Make own IGV where introns collapsed, output housekeeper for user
[ ] Cufflinks option to not do length norm
[ ] Pipeline test cases (paired end too)

Week ending June 9:
[ ] TCGA data
[ ] Ribosome profiling quantification with cufflinks needs GTF with only truncation is you want isoform abundance
[ ] UMI handling (add arguments and docs to trim and pipelines)
[ ] Meta-analysis test cases
  [ ] Make 3 nt light grid pattern to aid in viewing
  [ ] Weird that SRRXXXXX25 shows a drop at 5' end, look into this.
  [X] Complexity needs axis labels
  [ ] Give more notes as to how its quantifying (i.e. using longest transcript only?)
[ ] Work in Jon's tests
[ ] Strandedness protocol
[ ] Fix Travis and codecov and make sure all important steps are included in tests
[ ] Make sure aspects of XPRESSplot discussed are tested (PCA, volcano plots, normalizations)

Week ending June 16:
[ ] Test Docker instance
[ ] Figure out analysis output issues on HPC
  [ ] This is probably where the OOM errors on HPC are coming from. Works fine for smaller datasets, but not for bigger ones
[ ] Get test cases and docs feedback from people
[ ] Add note about truncator only editing exon files, thus must point htseq to exons
[ ] Use tximport to convert cufflinks to counts
[ ] Retry dataset with Cufflinks now that not masking and see how ISRIB stuff looks
[ ] Use a ribosome profiling specific DE package?
[ ] Incorporate Author feedback

Week ending June 23:
[ ] Send final manuscript to all authors for 1 week to get approval
[ ] Wrap up anything else for bioRXiv
[ ] Jupyter notebook examples
[ ] Walkthrough video
  - Add videos like example in this [README](https://github.com/manubot/manubot)
[ ] Benchmarking for Genome Biology?
[ ] Deconstruct diffs between Ingolia and XPRESSpipe
[ ] Submit to bioRXiv and Genome Biology

### After bioRXiv
[ ] Add instructions on how to add to the project
[ ] More tests!
[ ] Analyze other datasets using pipeline as backups


# Citation Notes
- DOI:10.1038/s41598-017-01617-3
  - Benchmarking of RNA-sequencing analysis workflows using whole transcriptome RT-qPCR expression data
  - All quantification methods about similar benchmarked to qPCR validation

# Other Notes
- NOTCH3 not annotated to be related to ISR from a quick google search
