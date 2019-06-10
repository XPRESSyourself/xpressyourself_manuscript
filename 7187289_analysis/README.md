# SLURM job 7187289
### Command used from XPRESSpipe https://github.com/XPRESSyourself/XPRESSpipe/commit/552b5351e9e1f78ceabf66426cdd855b2aa349be
```
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts_LCT.gtf -e isrib_riboprof -a CTGTAGGCACCATCAAT --method RPKM --sjdbOverhang 49 --quantification_method htseq
```
- New truncator that is fixed and parsed CDS
- Longest transcripts are Ensembl canonical transcripts
- Use htseq and riboseq module quantifies along CDS records (for both ribo and rna)
- No de-duplication for analysis
- Failed quality control execution due to OOM error 
