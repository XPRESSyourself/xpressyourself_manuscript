#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/u0690617/slurmjob-%j
#SBATCH --partition=kingspeak

### Running SRR1795425 and SRR1795426
# Ribo-seq Untreated Rep A Lane 1 & 2
# Possible issues
# -- rRNA masking
# -- de-duplication of reads in XPRESSpipe

### set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

FILES=/scratch/general/lustre/$USER/sra-files/walter_isrib_human_comparison/files
REF_INGOLIA=/scratch/general/lustre/$USER/sra-files/walter_isrib_human_comparison/references/ingolia
REF_STAR=/scratch/general/lustre/$USER/sra-files/walter_isrib_human_comparison/references/STAR
OUTPUT_INGOLIA=$SCRDIR/output/ingolia
OUTPUT_STAR=$SCRDIR/output/star_masking

mkdir $SCRDIR/input
mkdir $SCRDIR/output
mkdir $SCRDIR/output/ingolia
mkdir $SCRDIR/output/star_masking

cp $FILES/* $SCRDIR/input/.

cd $SCRDIR/.

### Prep references

# Make Ingolia refs


# Make STAR masking refs





### Pre-processing












### Run Ingolia
# Activate bowtie v1.0.0 & TopHat v2.0.7 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1606099)
conda activate ingolia_process ### Need to configure in KP to work


conda deactivate

### Run STAR


### Post-processing

# Run quantification on
