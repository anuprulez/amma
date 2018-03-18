source("https://bioconductor.org/biocLite.R")
suppressMessages(biocLite("DESeq2"))
suppressMessages(install.packages("pheatmap"))
suppressMessages(install.packages("gplots"))
suppressMessages(install.packages("UpSetR"))
suppressMessages(biocLite("org.Mm.eg.db"))
suppressMessages(install.packages("rentrez"))
suppressMessages(biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene"))

install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
suppressMessages(install.packages.2('devtools'))
suppressMessages(devtools::install_github('talgalili/dendextend', force=TRUE))

suppressMessages(biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")))
suppressMessages(install.packages("WGCNA"))
suppressMessages(biocLite("clusterProfiler"))
suppressMessages(install.packages('threejs'))