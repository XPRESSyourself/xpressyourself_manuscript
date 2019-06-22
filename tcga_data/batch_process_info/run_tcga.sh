#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/mammoth/serial/u0690617/slurmjob-%j
#SBATCH --partition=redwood-freecycle

PROJ=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab
REF=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab/references/human_reference_pe

# Process data
xpresspipe curateReference -o $REF -f $REF/genome_fastas --gtf $REF/transcripts.gtf --sjdbOverhang 75 -l

SCRDIR=/scratch/mammoth/serial/$USER/$SLURM_JOBID
mkdir -p $SCRDIR
cd $SCRDIR

FILES=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab/tcga/tcga_sample

mkdir $SCRDIR/input
mkdir $SCRDIR/output

cp $FILES/* $SCRDIR/input/.

cd $SCRDIR/.

# Process data
xpresspipe peRNAseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts_L.gtf -e tcga_validation -a None None --method FPKM --sjdbOverhang 75 --quantification_method htseq --two-pass

# Move output to project dir
zip -r $SCRDIR/tcga_validation.zip $SCRDIR/output
mv $SCRDIR/tcga_validation.zip $PROJ

rm -rf $SCRDIR



# Ran on 19_06_2019
# 17:26:44 Wall Clock   1-22:13:57 CPU Time   28 Cores    125Gn RAM
# Meta-analysis failed due to programmatic issues
