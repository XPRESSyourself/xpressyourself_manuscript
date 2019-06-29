#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/mammoth/serial/u0690617/slurmjob-%j
#SBATCH --partition=redwood-freecycle

# Process data
SCRDIR=/scratch/mammoth/serial/$USER/$SLURM_JOBID
mkdir -p $SCRDIR
cd $SCRDIR

PROJ=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab
REF=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab/references/human_reference_pe
FILES=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab/tcga/tcga_sample

mkdir $SCRDIR/input
mkdir $SCRDIR/output

cp $FILES/* $SCRDIR/input/.

cd $SCRDIR/.

xpresspipe curateReference -o $REF -f $REF/genome_fastas --gtf $REF/transcripts.gtf --sjdbOverhang 75

# Process data
xpresspipe peRNAseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts.gtf -e tcga_validation -a None None --sjdbOverhang 75 --quantification_method htseq --two-pass
