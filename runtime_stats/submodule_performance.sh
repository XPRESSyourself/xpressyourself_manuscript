#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/$USER/slurmjob-%j
#SBATCH --partition=kingspeak

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

SRA=/scratch/general/lustre/$USER/sra-files/walter_isrib_human
REF=/scratch/general/lustre/$USER/references/human_reference_se_v98

mkdir $SCRDIR/input
mkdir $SCRDIR/output

LOG=$SCRDIR/output/submodule_performance_stats.txt

cp $SRA/*.fastq $SCRDIR/input/

cd $SCRDIR/.


# Curate
echo "Start curate: $(date)" >> $LOG
/usr/bin/time -v xpresspipe makeReference -o $REF -f $REF/genome_fastas -g $REF/transcripts.gtf --sjdbOverhang 49
echo "End curate: $(date)" >> $LOG
echo "===================" >> $LOG


# Truncate GTF
echo "Start truncate: $(date)" >> $LOG
/usr/bin/time -v xpresspipe modifyGTF -g $REF/transcripts.gtf -t
echo "End truncate: $(date)" >> $LOG
echo "===================" >> $LOG


# Coding truncate GTF
echo "Start coding truncate: $(date)" >> $LOG
/usr/bin/time -v xpresspipe modifyGTF -g $REF/transcripts.gtf -p -t
echo "End coding truncate: $(date)" >> $LOG
echo "===================" >> $LOG


# Pipeline
echo "Start pipeline: $(date)" >> $LOG
/usr/bin/time -v xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/ -r $REF --gtf $REF/transcripts_CT.gtf -e test -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq
echo "End pipeline: $(date)" >> $LOG
echo "===================" >> $LOG


# Trim
echo "Start trim: $(date)" >> $LOG
/usr/bin/time -v xpresspipe trim -i $SCRDIR/input -o $SCRDIR/output -a CTGTAGGCACCATCAAT
echo "End trim: $(date)" >> $LOG
echo "===================" >> $LOG


# Align
echo "Start align: $(date)" >> $LOG
/usr/bin/time -v xpresspipe align -i $SCRDIR/output/trimmed_fastq -o $SCRDIR/output -r $REF -t SE --sjdbOverhang 49
echo "End align: $(date)" >> $LOG
echo "===================" >> $LOG


# abundance
echo "Start abundance: $(date)" >> $LOG
/usr/bin/time -v xpresspipe count -i $SCRDIR/output/alignments_coordinates -o $SCRDIR/output -g $REF/transcripts_CT.gtf --quantification_method cufflinks
echo "End abundance: $(date)" >> $LOG
echo "===================" >> $LOG

# Count
echo "Start count: $(date)" >> $LOG
/usr/bin/time -v xpresspipe count -i $SCRDIR/output/alignments_coordinates -o $SCRDIR/output -g $REF/transcripts_CT.gtf --feature_type CDS
echo "End count: $(date)" >> $LOG
echo "===================" >> $LOG

# DiffExpress
echo "Start diffex: $(date)" >> $LOG
cd /scratch/general/lustre/$USER/isrib_de/
/usr/bin/time -v bash run_de.sh
cd $SCRDIR/.
echo "End diffex: $(date)" >> $LOG
echo "===================" >> $LOG


# Read Distributions
echo "Start readDistribution: $(date)" >> $LOG
/usr/bin/time -v xpresspipe readDistribution -i $SCRDIR/output/trimmed_fastq -o $SCRDIR/output -t SE
echo "End readDistribution: $(date)" >> $LOG
echo "===================" >> $LOG


# Metagene
echo "Start Metagene: $(date)" >> $LOG
/usr/bin/time -v xpresspipe metagene -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf --feature_type CDS
echo "End Metagene: $(date)" >> $LOG
echo "===================" >> $LOG


# Gene Coverage
echo "Start Gene Coverage: $(date)" >> $LOG
/usr/bin/time -v xpresspipe geneCoverage -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf -n GAPDH
echo "End Gene Coverage: $(date)" >> $LOG
echo "===================" >> $LOG


# Periodicity
echo "Start Periodicity: $(date)" >> $LOG
/usr/bin/time -v xpresspipe periodicity -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf
echo "End Periodicity: $(date)" >> $LOG
echo "===================" >> $LOG


# Complexity
echo "Start Complexity: $(date)" >> $LOG
/usr/bin/time -v xpresspipe complexity -i $SCRDIR/output/alignments_coordinates -o $SCRDIR/output -g $REF/transcripts.gtf -t SE
echo "End Complexity: $(date)" >> $LOG
echo "===================" >> $LOG


# rRNA probe -- using pipeline output
echo "Start rRNA probe: $(date)" >> $LOG
/usr/bin/time -v xpresspipe rrnaProbe -i /scratch/general/lustre/u0690617/7559239/output/fastqc -o /scratch/general/lustre/u0690617/7559239/output
echo "End rRNA probe: $(date)" >> $LOG
echo "===================" >> $LOG
