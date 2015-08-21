# biomaRt_tutorial
A tutorial on using the biomaRt R package to retrieve data from online databases like Ensembl.

Simply follow along in biomaRt_demo.R and you will learn how to use the `biomaRt` package.

The GTF_converter.R script is an example of how to use `biomaRt` to do real work (converting gene ids in this case). The demo .gtf was created using the getGTF.sh shell script.

This tutorial is also available as a [webcast](https://plus.google.com/events/ca9v6dii91kapaqici7cmtqcgk4).

To install the required packages for this tutorial, enter the following commands in the R console:
```{r}
install.packages("stringr")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
```
