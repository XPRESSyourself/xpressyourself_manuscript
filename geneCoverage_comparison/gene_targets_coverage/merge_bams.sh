#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /scratch/general/lustre/$USER/slurmjob-%j
#SBATCH --partition=kingspeak

cd $SCRDIR/output/alignments_transcriptome/

# Merge transcriptome aligned bam files for geneCoverage

#samtools merge --threads 12 SRR1795409_1_Aligned.toTranscriptome.out.bam SRR1795410_1_Aligned.toTranscriptome.out.bam > SRR1795409_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795411_1_Aligned.toTranscriptome.out.bam SRR1795412_1_Aligned.toTranscriptome.out.bam > SRR1795411_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795413_1_Aligned.toTranscriptome.out.bam SRR1795414_1_Aligned.toTranscriptome.out.bam > SRR1795413_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795415_1_Aligned.toTranscriptome.out.bam SRR1795416_1_Aligned.toTranscriptome.out.bam > SRR1795415_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795417_1_Aligned.toTranscriptome.out.bam SRR1795418_1_Aligned.toTranscriptome.out.bam > SRR1795417_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795419_1_Aligned.toTranscriptome.out.bam SRR1795420_1_Aligned.toTranscriptome.out.bam > SRR1795419_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795421_1_Aligned.toTranscriptome.out.bam SRR1795422_1_Aligned.toTranscriptome.out.bam > SRR1795421_1_Merged.toTranscriptome.out.bam

#samtools merge --threads 12 SRR1795423_1_Aligned.toTranscriptome.out.bam SRR1795424_1_Aligned.toTranscriptome.out.bam > SRR1795423_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795425_1_Aligned.toTranscriptome.out.bam SRR1795426_1_Aligned.toTranscriptome.out.bam > SRR1795425_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795427_1_Aligned.toTranscriptome.out.bam SRR1795428_1_Aligned.toTranscriptome.out.bam > SRR1795427_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795429_1_Aligned.toTranscriptome.out.bam SRR1795430_1_Aligned.toTranscriptome.out.bam > SRR1795429_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795431_1_Aligned.toTranscriptome.out.bam SRR1795432_1_Aligned.toTranscriptome.out.bam > SRR1795431_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795433_1_Aligned.toTranscriptome.out.bam SRR1795434_1_Aligned.toTranscriptome.out.bam > SRR1795433_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795435_1_Aligned.toTranscriptome.out.bam SRR1795436_1_Aligned.toTranscriptome.out.bam > SRR1795435_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795437_1_Aligned.toTranscriptome.out.bam SRR1795438_1_Aligned.toTranscriptome.out.bam > SRR1795437_1_Merged.toTranscriptome.out.bam

samtools merge --threads 12 SRR1795439_1_Aligned.toTranscriptome.out.bam SRR1795440_1_Aligned.toTranscriptome.out.bam > SRR1795439_1_Merged.toTranscriptome.out.bam

mkdir merged_bams
mv *_Merged*bam merged_bams/


# POMGNT1, RPL27, TKT, HSPA8, NDUFA11
for X in POMGNT1 RPL27 TKT HSPA8 NDUFA11;
do
  GENE=${X};
  mkdir $GENE;
  xpresspipe geneCoverage -i merged_bams/ -o $GENE/ -g transcripts.gtf -n $GENE \
  --samples SRR1795425_1_Merged.toTranscriptome.out.bam SRR1795427_1_Merged.toTranscriptome.out.bam SRR1795433_1_Merged.toTranscriptome.out.bam SRR1795435_1_Merged.toTranscriptome.out.bam SRR1795437_1_Merged.toTranscriptome.out.bam SRR1795439_1_Merged.toTranscriptome.out.bam SRR1795429_1_Merged.toTranscriptome.out.bam SRR1795431_1_Merged.toTranscriptome.out.bam \
  --sample_names ribo_untr_a ribo_untr_b ribo_tm_a ribo_tm_b ribo_tmisrib_a ribo_tmisrib_b ribo_isrib_a ribo_isrib_b;
done
