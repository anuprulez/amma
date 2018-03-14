source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("pheatmap")
install.packages("gplots")
install.packages("UpSetR")
biocLite("org.Mm.eg.db")
install.packages("rentrez")
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")

install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
install.packages.2('devtools')
devtools::install_github('talgalili/dendextend', force=TRUE)

biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")

biocLite("clusterProfiler")

install.packages('threejs')