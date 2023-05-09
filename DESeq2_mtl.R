
#####################
# 	Load libraries 	#
#####################
library("DESeq2")

#####################
# 	Load data 	#
#####################
exprmatrix<- read.table("ASE_pop_deseq2_countMatrix_genes_14Jan2020_M", header=FALSE)
exprmatrix[is.na(exprmatrix)] <- 0

sampletable<- read.table("ASE_pop_deseq2_sampletable_genes_7Jan2020_M_uniq_nohead", header=FALSE)
samps<- data.frame(sampletable)
colnames(samps)<- c("individual", "allele")

#####################
# 	Create matrix 	#
#####################
bigmat<- data.matrix(exprmatrix[2:ncol(exprmatrix)])
mode(bigmat) <- "integer"
colnames(bigmat)<- samps$individual
gene<- exprmatrix[1]

#####################
# 	Run DESeq2 	#
#####################

design <- model.matrix(~0+samps$individual+ samps$allele)

dds <- DESeqDataSetFromMatrix(countData = bigmat, colData = samps, design = design)

sizeFactors(dds) <- rep(1, length(samps$individual))
dds <- DESeq(dds)

res<-results(dds, name="samps.alleleDer")
res<- data.frame(res)
res2<-merge(gene, res, by=c("row.names"), all=TRUE)

write.table(res2, "DeSeq2_Results_alleleDer_M_14Jan2020", sep="\t", row.names=FALSE, col.names=TRUE)
