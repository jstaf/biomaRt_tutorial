#!/bin/bash

# grab our data from the ensembl ftp site
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-26/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.26.gff3.gz
gunzip *.gz

# you can see whats in the file with the following command
# awk -F "\t" '{print $3}' Drosophila_melanogaster.BDGP6.26.gff3 | sort | uniq

# It's normally somewhat hard to convert a file to GTF from GFF3, but the gffread program from cufflinks (http://cole-trapnell-lab.github.io/cufflinks/) can do it for us.
gffread -E Drosophila_melanogaster.BDGP6.26.gff3 -T -o dmel_geneset-r6.26.gtf

# if you re-run the awk bit from above, youll notice that now we only have pseudogenes, exons, and CDS

# There's unfortunately some residual bad stuff in the gene_id and transcript_id columns that we'll need to get rid of.
# sed can do a find and replace on the weird bits
sed -i 's/transcript://g' dmel_geneset-r6.26.gtf
sed -i 's/gene://g' dmel_geneset-r6.26.gtf

# All done! We have our starting GTF.
