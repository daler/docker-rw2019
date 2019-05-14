# Demonstration of running: create a DESeq2 example dataset and save a PDF of
# the PCA of those data
library(DESeq2)
dds <- makeExampleDESeqDataSet(betaSD=1)
rld <- rlog(dds)
pdf(file='pca.pdf')
plotPCA(rld)
dev.off()
sessionInfo()
