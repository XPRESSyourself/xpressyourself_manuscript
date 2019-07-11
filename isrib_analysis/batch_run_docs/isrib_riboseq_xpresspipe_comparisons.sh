#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/uNID/slurmjob-%j
#SBATCH --partition=kingspeak

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

SRA=/scratch/general/lustre/$USER/sra-files/walter_isrib_human
REF2=/scratch/general/lustre/uNID/references/human_reference_se

mkdir $SCRDIR/input
mkdir $SCRDIR/output

cp $SRA/*.fastq $SCRDIR/input/.

cd $SCRDIR/.

xpresspipe curateReference -o $REF2 -f $REF2/genome_fastas -g $REF2/transcripts.gtf -l -p -t --sjdbOverhang 49
xpresspipe modifyGTF -g $REF2/transcripts.gtf -p -t

"""slurmjob-7322864"""
mkdir $SCRDIR/output/v96_longest_truncated
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v96_longest_truncated -r $REF2 --gtf $REF2/transcripts_LCT.gtf -e isrib_comp_v96_longest_truncated -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq

"""slurmjob-7322865"""
mkdir $SCRDIR/output/v96_truncated
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v96_truncated -r $REF2 --gtf $REF2/transcripts_CT.gtf -e isrib_comp_v96_truncated -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq

"""slurmjob-7322866"""
mkdir $SCRDIR/output/v96_normal
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v96_normal -r $REF2 --gtf $REF2/transcripts.gtf -e isrib_comp_v96_normal -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq



"""slurmjob-7325061"""
REF=/scratch/general/lustre/uNID/references/human_reference_se_v72

mkdir $SCRDIR/input
mkdir $SCRDIR/output

cp $SRA/*.fastq $SCRDIR/input/.

cd $SCRDIR/.

xpresspipe curateReference -o $REF -f $REF/genome_fastas -g $REF/transcripts.gtf -l -p -t --sjdbOverhang 49
xpresspipe modifyGTF -g $REF/transcripts.gtf -p -t

#run script --- add directory of raw data
#Using Qiagen Illumina 3' adaptor, but shouldn't need adaptor trimming anyways, just a placeholder for now
mkdir $SCRDIR/output/v72_longest_truncated
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v72_longest_truncated -r $REF --gtf $REF/transcripts_LCT.gtf -e isrib_comp_v72_longest_truncated -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq


mkdir $SCRDIR/output/v72_truncated
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v72_truncated -r $REF --gtf $REF/transcripts_CT.gtf -e isrib_comp_v72_truncated -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq


mkdir $SCRDIR/output/v72_normal
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v72_normal -r $REF --gtf $REF/transcripts.gtf -e isrib_comp_v72_normal -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq
