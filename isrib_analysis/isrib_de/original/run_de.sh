xpresspipe diffxpress -i tm_counts.txt -s tm_deseq.txt --design Type+Condition+Type:Condition
xpresspipe diffxpress -i tmisrib_counts.txt -s tmisrib_deseq.txt --design Type+Condition+Type:Condition
xpresspipe diffxpress -i isrib_counts.txt -s isrib_deseq.txt --design Type+Condition+Type:Condition

xpresspipe diffxpress -i tm_ribo_counts.txt -s tm_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i tm_rna_counts.txt -s tm_rna_deseq.txt --design Condition
xpresspipe diffxpress -i tmisrib_ribo_counts.txt -s tmisrib_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i tmisrib_rna_counts.txt -s tmisrib_rna_deseq.txt --design Condition
xpresspipe diffxpress -i isrib_ribo_counts.txt -s isrib_ribo_deseq.txt --design Condition
xpresspipe diffxpress -i isrib_rna_counts.txt -s isrib_rna_deseq.txt --design Condition
