#!/bin/bash

# Analyze differences in two recent releases of Homo Sapiens GTF file to highlight rapid evolution of annotations and how might effect TCGA comparisons between XPRESSpipe and their processing pipeline

# Get GTFs (21 Jun 2019, 11:40AM)
# Last Modified 11/24/18, 5:00:00 PM
curl -O ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
# Last Modified 3/13/19, 6:03:00 AM
curl -O ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
gzip -d *gtf.gz

# Get line stats
cat cat Homo_sapiens.GRCh38.95.gtf | wc -l
# 2737564
cat Homo_sapiens.GRCh38.96.gtf | wc -l
# 2769823
expr 2769823 - 2737564
# 32259
