xpresspipe diffxpress -i tm_counts.tsv -s tm_deseq.txt --design Type+Condition+Type:Condition
xpresspipe diffxpress -i tmisrib_counts.tsv -s tmisrib_deseq.txt --design Type+Condition+Type:Condition
xpresspipe diffxpress -i isrib_counts.tsv -s isrib_deseq.txt --design Type+Condition+Type:Condition

xpresspipe diffxpress -i tm_ribo_counts.tsv -s tm_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i tm_rna_counts.tsv -s tm_rna_deseq.txt --design Condition
xpresspipe diffxpress -i tmisrib_ribo_counts.tsv -s tmisrib_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i tmisrib_rna_counts.tsv -s tmisrib_rna_deseq.txt --design Condition
xpresspipe diffxpress -i isrib_ribo_counts.tsv -s isrib_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i isrib_rna_counts.tsv -s isrib_rna_deseq.txt --design Condition
