#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=11
#SBATCH -o /scratch/mammoth/serial/u0690617/slurmjob-%j
#SBATCH --partition=redwood-freecycle

PROJ=/uufs/chpc.utah.edu/common/HIPAA/proj_rutter_lab
DIR=$PROJ/tcga_comp
REF_79=$PROJ/tcga_comp/ref_79
GTF_79=$PROJ/tcga_comp/ref_79/transcripts.gtf
GTF_79L=$PROJ/tcga_comp/ref_79/transcripts_L.gtf
REF_96=$PROJ/tcga_comp/ref_96
GTF_96=$PROJ/tcga_comp/ref_96/transcripts.gtf
GTF_96L=$PROJ/tcga_comp/ref_96/transcripts_L.gtf
FILE_DIR=$PROJ/tcga_comp/files


# Prep references
xpresspipe curateReference -o $REF_79 -f $REF_79 -g $GTF_79 -l --sjdbOverhang 100
xpresspipe curateReference -o $REF_96 -f $REF_96 -g $GTF_96 -l --sjdbOverhang 100


# Trim strict
# Index v79, GTF and fasta
# Unique only -- strict
# Quant v79 whole GTF
mkdir $DIR/v79_strictTrim_normalGTF_uniqueOnly;
mkdir $DIR/v79_strictTrim_normalGTF_multimappers
mkdir $DIR/v79_strictTrim_longestGTF_uniqueonly
mkdir $DIR/v79_strictTrim_longestGTF_multimappers
mkdir $DIR/v96_strictTrim_normalGTF_uniqueOnly
mkdir $DIR/v96_strictTrim_normalGTF_multimappers
mkdir $DIR/v96_strictTrim_longestGTF_uniqueonly
mkdir $DIR/v96_strictTrim_longestGTF_multimappers
mkdir $DIR/v79_looseTrim_normalGTF_uniqueOnly
mkdir $DIR/v79_looseTrim_normalGTF_multimappers
mkdir $DIR/v79_looseTrim_longestGTF_uniqueonly
mkdir $DIR/v79_looseTrim_longestGTF_multimappers
mkdir $DIR/v96_looseTrim_normalGTF_uniqueOnly
mkdir $DIR/v96_looseTrim_normalGTF_multimappers
mkdir $DIR/v96_looseTrim_longestGTF_uniqueonly
mkdir $DIR/v96_looseTrim_longestGTF_multimappers


srun -n 1 -c 16 xpresspipe peRNAseq \
    -i $FILE_DIR \
    -o $DIR/v79_strictTrim_normalGTF_uniqueOnly \
    -e v79_strictTrim_normalGTF_uniqueOnly \
    -r $REF_79 \
    -g $GTF_79 \
    --sjdbOverhang 100 \
    -a None None \
    -c htseq \
    --two-pass & \



srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_strictTrim_normalGTF_multimappers \
  -e v79_strictTrim_normalGTF_multimappers \
  -r $REF_79 \
  -g $GTF_79 \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v79, GTF and fasta
# Multimappers allowed
# Quant v79 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_strictTrim_longestGTF_uniqueonly \
  -e v79_strictTrim_longestGTF_uniqueonly \
  -r $REF_79 \
  -g $GTF_79L \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v79, GTF and fasta
# Multimappers allowed
# Quant v79 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_strictTrim_longestGTF_multimappers \
  -e v79_strictTrim_longestGTF_multimappers \
  -r $REF_79 \
  -g $GTF_79L \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v96, GTF and fasta
# Unique only -- strict
# Quant v96 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_strictTrim_normalGTF_uniqueOnly \
  -e v96_strictTrim_normalGTF_uniqueOnly \
  -r $REF_96 \
  -g $GTF_96 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v96, GTF and fasta
# Unique only -- strict
# Quant v96 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_strictTrim_normalGTF_multimappers \
  -e v96_strictTrim_normalGTF_multimappers \
  -r $REF_96 \
  -g $GTF_96 \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v96, GTF and fasta
# Multimappers allowed
# Quant v96 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_strictTrim_longestGTF_uniqueonly \
  -e v96_strictTrim_longestGTF_uniqueonly \
  -r $REF_96 \
  -g $GTF_96L \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim strict
# Index v96, GTF and fasta
# Multimappers allowed
# Quant v96 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_strictTrim_longestGTF_multimappers \
  -e v96_strictTrim_longestGTF_multimappers \
  -r $REF_96 \
  -g $GTF_96L \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \






# Trim loose
# Index v79, GTF and fasta
# Unique only -- strict
# Quant v79 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_looseTrim_normalGTF_uniqueOnly \
  -e v79_looseTrim_normalGTF_uniqueOnly \
  -r $REF_79 \
  -g $GTF_79 \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v79, GTF and fasta
# Unique only -- strict
# Quant v79 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_looseTrim_normalGTF_multimappers \
  -e v79_looseTrim_normalGTF_multimappers \
  -r $REF_79 \
  -g $GTF_79 \
  --allow_multimappers \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v79, GTF and fasta
# Multimappers allowed
# Quant v79 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_looseTrim_longestGTF_uniqueonly \
  -e v79_looseTrim_longestGTF_uniqueonly \
  -r $REF_79 \
  -g $GTF_79L \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v79, GTF and fasta
# Multimappers allowed
# Quant v79 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v79_looseTrim_longestGTF_multimappers \
  -e v79_looseTrim_longestGTF_multimappers \
  -r $REF_79 \
  -g $GTF_79L \
  -q 10 \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v96, GTF and fasta
# Unique only -- strict
# Quant v96 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_looseTrim_normalGTF_uniqueOnly \
  -e v96_looseTrim_normalGTF_uniqueOnly \
  -r $REF_96 \
  -g $GTF_96 \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v96, GTF and fasta
# Unique only -- strict
# Quant v96 whole GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_looseTrim_normalGTF_multimappers \
  -e v96_looseTrim_normalGTF_multimappers \
  -r $REF_96 \
  -g $GTF_96 \
  --allow_multimappers \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v96, GTF and fasta
# Multimappers allowed
# Quant v96 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_looseTrim_longestGTF_uniqueonly \
  -e v96_looseTrim_longestGTF_uniqueonly \
  -r $REF_96 \
  -g $GTF_96L \
  -q 10 \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \


# Trim loose
# Index v96, GTF and fasta
# Multimappers allowed
# Quant v96 longest GTF

srun -n 1 -c 16 xpresspipe peRNAseq \
  -i $FILE_DIR \
  -o $DIR/v96_looseTrim_longestGTF_multimappers \
  -e v96_looseTrim_longestGTF_multimappers \
  -r $REF_96 \
  -g $GTF_96L \
  -q 10 \
  --allow_multimappers \
  --sjdbOverhang 100 \
  -a None None \
  -c htseq \
  --two-pass & \

wait




mv $DIR/v79_strictTrim_normalGTF_uniqueOnly/counts/* $DIR
mv $DIR/v79_strictTrim_normalGTF_multimappers/counts/* $DIR
mv $DIR/v79_strictTrim_longestGTF_uniqueonly/counts/* $DIR
mv $DIR/v79_strictTrim_longestGTF_multimappers/counts/* $DIR
mv $DIR/v96_strictTrim_normalGTF_uniqueOnly/counts/* $DIR
mv $DIR/v96_strictTrim_normalGTF_multimappers/counts/* $DIR
mv $DIR/v96_strictTrim_longestGTF_uniqueonly/counts/* $DIR
mv $DIR/v96_strictTrim_longestGTF_multimappers/counts/* $DIR
mv $DIR/v79_looseTrim_normalGTF_uniqueOnly/counts/* $DIR
mv $DIR/v79_looseTrim_normalGTF_multimappers/counts/* $DIR
mv $DIR/v79_looseTrim_longestGTF_uniqueonly/counts/* $DIR
mv $DIR/v79_looseTrim_longestGTF_multimappers/counts/* $DIR
mv $DIR/v96_looseTrim_normalGTF_uniqueOnly/counts/* $DIR
mv $DIR/v96_looseTrim_normalGTF_multimappers/counts/* $DIR
mv $DIR/v96_looseTrim_longestGTF_uniqueonly/counts/* $DIR
mv $DIR/v96_looseTrim_longestGTF_multimappers/counts/* $DIR
