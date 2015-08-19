# Let's convert gene ids in a GTF file from FBgn#s to human intelligible names
# Jeff Stafford

# What is a GTF file? It's a file format commonly used to describe the locations
# of genetic features in a genome (its often used in sequencing and heavy duty 
# bioinformatics applications). These files get huge, and can contain lots of 
# annotations, however the format is extremely strict, so using them always 
# happens the same way. It also happens to be a good example of where biomaRt is
# useful.

# GTF Format Specification: 

# Files are tab delimited, with each column having a specific type of data
# Column 1 = Sequence name (contig # or chromosome name)
# Column 2 = Source - where did this data come from?
# Column 3 = What type of feature is this? Is it a gene/CDS/exon/stop codon/etc?
# Column 4 = Sequence start (beginning of a chromosome is "1")
# Column 5 = Sequence end
# Column 6 = Score - confidence of the annotation's accuracy. Can just be a "." 
# Column 7 = Strand - Which DNA strand is it on? (+ or -)
# Column 8 = Frame - What frame is the sequence in? (0, 1, or 2. "." = NA)
# Column 9 = Metadata... this one's tricky. Metadata is actually another set of
# columns that are delimited by ";"s. It generally looks like this:
# gene_id "geneName"; transcript_id "transcriptName"; moreData "Value"; 

# Anyhow, how do we load that in R? Seems tricky, right?

# Our data for today is every exon, coding sequence, and pseudogene in the 
# Drosophila melanogaster genome... this file was obtained as a GFF3 from
# Ensembl Metazoa using the getGTF.sh script in this repo.

# load a .gtf as a dataframe... 
gtf <- read.delim("dmel_geneset_abridged-r6.26.gtf",
                  header = FALSE, as.is = TRUE, quote = "")
# we have to use the 'quote = ""' bit or R will get rid of our quotation marks
# automatically, which sucks.
str(gtf)

# We are interested only in grabbing the gene ID from the 9th column
# metadata is contained in column 9
metadata <- gtf[, 9]
head(metadata)

# split the data by \"
library(stringr)
metadata <- str_split(metadata, pattern = "\"")
# put metadata back into matrix of strings instead of list of strings
metadata <- do.call(rbind, metadata)

# now get unique gene id's and retrieve the actual gene names from ensembl using biomaRt
library(biomaRt)
uniq_genes <- unique(metadata[, 4])
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl") # connect to a mart
martData <- getBM(attributes = c("ensembl_gene_id", # fbgn#s work as ensembl ids
                                 "flybasename_gene"), # actual gene names
                  filters = "ensembl_gene_id", # filters = type of values we are searching by
                  values = uniq_genes, # values to search with
                  mart = ensembl) # mart we are searching
head(martData) # yay we have our gene names

# match our data to the retrieved data and replace...
# protip with match: 1st arg = what we are matching to, 2nd arg = what we want to match
idx <- match(metadata[, 4], martData$ensembl_gene_id)
replacement <- martData$flybasename_gene[idx]
# notice the NAs, we don't want to replace those ones where we couldnt grab data
head(replacement)
notNA <- which(!is.na(replacement))
metadata[notNA, 4] <- replacement[notNA]
# now we can see that the data is fixed
head(metadata)

# okay now put everything back together and export, remember that we need those quotation marks!
quotes <- rep("\"", length(metadata[, 1]))
metadata <- paste(metadata[, 1], quotes, metadata[, 2], quotes, metadata[, 3],
              quotes, metadata[, 4], quotes, metadata[, 5],
              sep = "")
# check our work
head(metadata)

# Looks good! Let's export it!
gtf[, 9] <- metadata
write.table(gtf,
            file = "dmel_geneset_relabel-r6.26.gtf",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

