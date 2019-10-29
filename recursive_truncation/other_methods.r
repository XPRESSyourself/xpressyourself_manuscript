# A brief look at some alternative methods for ignoring ends of CDSs for Quantification

library(GenomicFeatures)
library(GenomicAlignments)
library(data.table)

gtf_file <- "Saccharomyces_cerevisiae.R64-1-1.97.gtf"
bam_file <- "SRR1795435_1_Merged.toTranscriptome.out.bam"

# Modifying GTF using R
# Could add 45 to l_utr5 and subtract 45 from l_cds, but won't account for edge cases
# Also, requires proficiency in R, no way to export, etc
setwd("~/Desktop/trunc_test/")
gtf <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file, format = 'gtf')

gene <- suppressWarnings(GenomicFeatures::transcriptsBy(gtf, "gene"))
exon <- suppressWarnings(GenomicFeatures::exonsBy(gtf, by = "tx",use.names=T))
utr5 <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(gtf,use.names=T))
cds <- suppressWarnings(GenomicFeatures::cdsBy(gtf, by = "tx", use.names=T))
utr3 <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(gtf,use.names=T))

# Convert features to dataframes
gene <- as.data.table(gene)
exon <- as.data.table(exon[unique(names(exon))])
utr5 <- as.data.table(utr5[unique(names(utr5))])
cds <- as.data.table(cds[unique(names(cds))])
utr3 <-as.data.table(utr3[unique(names(utr3))])

# Parse out relevant information, with transcript_id being the key
gene_name <- gene[, c('tx_name','group_name')]
names(gene_name) <- c("transcript", "gene")
l_transcript <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]

# Merge records
print('Collecting flat reference information...')

merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
flat_gtf  <-  Reduce(merge_allx, list(gene_name, l_transcript, l_utr5, l_cds, l_utr3))
flat_gtf[is.na(flat_gtf)] <- 0

rows <- dim(flat_gtf)[1]
print(paste(rows, 'transcript records compiled...', sep=' '))

head(flat_gtf)

# Modifying transcriptome alignments in R
# Requires transcriptome alignment, less flexible for downstream purposes
# Assumes the user is proficient in R, can account for edge cases, and can determine appropriately when to trim in this context
# No information on where CDS is here
bam <- as.data.table(GenomicAlignments::readGAlignments(bam_file))
head(bam)
