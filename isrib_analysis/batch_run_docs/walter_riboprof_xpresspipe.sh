#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/u0690617/slurmjob-%j
#SBATCH --partition=kingspeak

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

SRA=/scratch/general/lustre/$USER/sra-files/walter_isrib_human
REF=/scratch/general/lustre/u0690617/references/human_reference_se

mkdir $SCRDIR/input
mkdir $SCRDIR/output

cp $SRA/*.fastq $SCRDIR/input/.

cd $SCRDIR/.

#run script --- add directory of raw data
#Using Qiagen Illumina 3' adaptor, but shouldn't need adaptor trimming anyways, just a placeholder for now
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts_LCT.gtf -e isrib_riboprof -a CTGTAGGCACCATCAAT --method RPM --sjdbOverhang 49 --quantification_method htseq
