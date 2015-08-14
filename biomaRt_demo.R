# So... what is biomaRt?
library(biomaRt)

# While many R packages focus on working with data that you have, biomaRt can 
# grab you the data that you don't have. Specifically, it can connect to online
# bioinformatics databases and grab data from them.

# NOTE: YOU NEED TO BE CONNECTED TO THE INTERNET FOR BIOMART TO WORK.

# Which databases are available?
marts <- listMarts()
head(marts)

# There are a lot of different databases out there... some, like Ensembl, have 
# multiple versions for different stuff. Others, like TAIR, focus on a 
# particular organism (like Arabidopsis). Note that occasionally the mart names
# change from time to time as stuff gets updated.

# We're going to use ensembl genes.
ensembl <- useMart("ensembl")

# What datasets does ensembl have in it? 
datasets <- listDatasets(ensembl)
head(datasets)

# We're going to use the Drosophila melanogaster stuff, because that's what I
# work on.
ensembl <- useMart("ensembl", "dmelanogaster_gene_ensembl")

# What data is available for us?
attributes <- listAttributes(ensembl)
# There's a lot... over 1128 different things we can retrieve!
dim(attributes)
head(attributes)

# So, how do we actually get things from ensembl in this case?

# Here's a random list of Drosophila genes. I've only heard of one of them
# before, and have no idea what most of the others do. Lets find out!
geneList <- read.csv("gene_list.csv")
geneList <- geneList[,2]

# getBM() is the workhorse for biomaRt and is used to retrieve data.

# attributes are what we want to retrieve... they're kind of limitless
martData <- getBM(attributes = c("ensembl_gene_id", # FBgn#, highly specific identifier, changes often
                                 "flybasename_gene", # the normal gene name
                                 "go_id", # go term accession #
                                 "name_1006", # go term name
                                 "transcript_count", # number of transcripts
                                 "transmembrane_domain"), # transmembrane domain type (why not?)
                  filters = "flybasename_gene", # filters = type of values we are searching by
                  values = geneList, # values to search with
                  mart = ensembl) # mart we are searching

# Note: you can't use all attributes as filters... to see which ones you can use
# you can use listFilters().

# got all of the data back, but there are multiple entries for items when we
# retrieved multiple values for an gene (like for GO terms)
str(martData) 
length(unique(martData$flybasename_gene)) # we found data for all 20 genes

# A lot to take in, isn't it?

# It looks like Galphaq is a GTPase (GO:0003924- GTPase activity). What if we
# randomly REALLY wanted to know what other genes were involved in that? Is
# there a way to retrieve every other gene known to be a GTPase?
GTPases <- getBM(attributes = c("ensembl_gene_id", "flybasename_gene", "go_id"),
                 filters = "go_id",
                 values = "GO:0003924",
                 mart = ensembl)
# This should be every annotated GTPase in the Drosophila genome
head(GTPases)

# What if we wanted the transcript nucleotide sequence for every transcript from
# our list of genes?
seqs <- getSequence(id = geneList, # same as values from getBM()
                    type = "flybasename_gene", # same as filters
                    seqType = "transcript_exon", 
                    mart = ensembl)
# Notice that once again, we've retrieved multiple entries for each gene, how do
# we know what sequence corresponds to what transcript?
str(seqs)

# What if we wanted to be smarter about it, and link every transcript to it's
# transcript name?

# First, lets get all of the different transcripts.
transcripts <- getBM(attributes = c("ensembl_gene_id", 
                                    "flybasename_gene", 
                                    "flybase_transcript_id"),
                     filters = "flybasename_gene",
                     values = geneList,
                     mart = ensembl)
# Note that biomaRt sorts output alphabetically.
str(transcripts) 
# Now lets get the sequences. 
seqs <- getSequence(id = transcripts$flybase_transcript_id, 
                    type = "flybase_transcript_id", 
                    seqType = "peptide",
                    mart = ensembl)
# Lets match the data in "seqs" to the gene names in "transcripts"
# match() is tricky
# first arg = what we are matching to
# second arg = what we want to match
index <- match(transcripts$flybase_transcript_id, seqs$flybase_transcript_id)

# we can use the index to match things up correctly
transcripts$peptide_sequence <- seqs$peptide[index]
# Now lets sort the list of transcripts and genes alphabetically

# Beautiful. Now lets write these sequences to disk as fasta files just for kicks.

# Make a folder
dir.create("fasta_seqs")

# I'm going to write one file with all of the peptide sequences for each gene. 
# It looks like the biomaRt devs gave us a nice exportFASTA function. However, 
# it only works with the raw output of getSequence(), so we are going to make a
# function to put all transcripts in a file, one per gene.

makeFASTA <- function(gene) {
  filename <- paste("fasta_seqs/", gene, ".fa", sep = "")
  trans <- subset(transcripts$flybase_transcript_id, transcripts$flybasename_gene == gene)
  sequences <- subset(seqs, seqs$flybase_transcript_id %in% trans)
  exportFASTA(sequences, filename)
  return(paste("Exporting", gene))
}
# now lets run this across all genes
lapply(unique(transcripts$flybasename_gene), makeFASTA)
