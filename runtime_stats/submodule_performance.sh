#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/$USER/slurmjob-%j
#SBATCH --partition=kingspeak

#set up the temporary directory
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

SRA=/scratch/general/lustre/$USER/sra-files/walter_isrib_human
REF=/scratch/general/lustre/$USER/references/human_reference_se

mkdir $SCRDIR/input
mkdir $SCRDIR/output

LOG=$SCRDIR/output/submodule_performance_stats.txt

cp $SRA/*.fastq $SCRDIR/input/

cd $SCRDIR/.


# mkdir $SCRDIR/output/v76_longest_truncated
# xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output/v76_longest_truncated -r $REF --gtf $REF/transcripts_LCT.gtf -e # isrib_comp_v76_longest_truncated -a CTGTAGGCACCATCAAT --sjdbOverhang 49 --quantification_method htseq


# Curate
echo "Start curate: $(date)" >> $LOG
xpresspipe makeReference -o $REF -f $REF/genome_fastas -g $REF/transcripts.gtf --sjdbOverhang 49
echo "End curate: $(date)" >> $LOG
echo "===================" >> $LOG


# Truncate GTF
echo "Start truncate: $(date)" >> $LOG
xpresspipe modifyGTF -g $REF/transcripts.gtf -t
echo "End truncate: $(date)" >> $LOG
echo "===================" >> $LOG


# Trim
echo "Start trim: $(date)" >> $LOG
xpresspipe trim -i $SCRDIR/input -o $SCRDIR/output -a CTGTAGGCACCATCAAT
echo "End trim: $(date)" >> $LOG
echo "===================" >> $LOG


# Align
echo "Start align: $(date)" >> $LOG
xpresspipe align -i $SCRDIR/output/trimmed_fastq -o $SCRDIR/output -r $REF -t SE --sjdbOverhang 49
echo "End align: $(date)" >> $LOG
echo "===================" >> $LOG


# Count
echo "Start count: $(date)" >> $LOG
xpresspipe count -i $SCRDIR/output/alignments_coordinates -o $SCRDIR/output -g $REF/transcripts_T.gtf --feature_type CDS
echo "End count: $(date)" >> $LOG
echo "===================" >> $LOG


# DiffExpress
echo "Start diffex: $(date)" >> $LOG
cd /scratch/general/lustre/$USER/isrib_de/
bash run_de.sh
cd $SCRDIR/.
echo "End diffex: $(date)" >> $LOG
echo "===================" >> $LOG


# Read Distributions
echo "Start readDistribution: $(date)" >> $LOG
xpresspipe readDistribution -i $SCRDIR/output/trimmed_fastq -o $SCRDIR/output -t SE
echo "End readDistribution: $(date)" >> $LOG
echo "===================" >> $LOG


# Metagene
echo "Start Metagene: $(date)" >> $LOG
xpresspipe metagene -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf --feature_type CDS
echo "End Metagene: $(date)" >> $LOG
echo "===================" >> $LOG


# Gene Coverage
echo "Start Gene Coverage: $(date)" >> $LOG
xpresspipe geneCoverage -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf -n GAPDH
echo "End Gene Coverage: $(date)" >> $LOG
echo "===================" >> $LOG


# Periodicity
echo "Start Periodicity: $(date)" >> $LOG
xpresspipe periodicity -i $SCRDIR/output/alignments_transcriptome -o $SCRDIR/output -g $REF/transcripts.gtf
echo "End Periodicity: $(date)" >> $LOG
echo "===================" >> $LOG


# Complexity
echo "Start Complexity: $(date)" >> $LOG
xpresspipe complexity -i $SCRDIR/output/alignments_coordinates -o $SCRDIR/output -g $REF/transcripts.gtf -t SE
echo "End Complexity: $(date)" >> $LOG
echo "===================" >> $LOG


# rRNA probe
echo "Start rRNA probe: $(date)" >> $LOG
mkdir $SCRDIR/output/fastqc
for X in $SCRDIR/input/trimmed_fastq/*.fastq; do fastqc -q $SCRDIR/input/trimmed_fastq/${X} -o $SCRDIR/output/fastqc; done
xpresspipe rrnaProbe -i $SCRDIR/output/fastqc -o $SCRDIR/output
echo "End rRNA probe: $(date)" >> $LOG
echo "===================" >> $LOG
