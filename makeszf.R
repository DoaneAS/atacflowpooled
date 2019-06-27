

library(DESeq2)
library(tidyverse)

ff =  Sys.glob("*txt.gz")
ff
nb = read_tsv("nbcounts.counts.txt.gz", col_names="nb")


gcb = read_tsv("gcbcounts.counts.txt.gz", col_names="gcb")

ly7wt = read_tsv("LY7WT.counts.txt.gz", col_names="ly7wt")

ly7ko = read_tsv("LY7_OCT2KO.counts.txt.gz", col_names="ly7OCT2KO")

cts = cbind.data.frame(nb, gcb, ly7wt, ly7ko)



load("/athena/elementolab/scratch/asd2007/projectshg38/analysis/bcell/cellAtlasB.Rdata")

biocLite("DESeq2")

ctss =  cbind(assay(cellAtlasB_RSE), cts)


colnames(ctss)

pd

library(SummarizedExperiment)



atlas = rowRanges(cellAtlasB_RSE)

pd = colData(cellAtlasB_RSE)

pdd = data.frame(sampleID=colnames(ctss), cellType=c(as.character(pd$cellTypeg), "NB", "GCB", "LY7_WT", "LY7_OCT2KO"))

str(ctss)

ctss = data.matrix(ctss)

sef =  SummarizedExperiment(assays=list(counts=ctss), colData=pdd)

rowRanges(sef) = atlas

pd
atlas

atlas = reduce(atlas)

#pdd  = rbind(pd, pd[3,], pd[1,])

#rownames(pdd)[26] <- "gcb"

#rownames(pdd)[25] <- "nb"

d = DESeqDataSet(sef, design= ~cellType)

dd =estimateSizeFactors(d)
norms =  1/sizeFactors(dd)

nf = data.frame(sampleID =  names(norms), sf=as.vector(norms))


write_tsv(nf, "sizefactors.txt")

sef
rowRanges(sef)

ddsSE <- DESeqDataSet(sef, design = ~ cellType )

load("sef.Rdata")
save(sef, file="sef.Rdata")

kk
