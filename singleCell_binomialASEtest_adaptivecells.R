#####################
# 	Load libraries 	#
#####################
library(dplyr)
library(stringr)
library(ggplot2)

#####################
# 	Load data 	#
#####################
cell<- read.table("/Users/mharwood/Documents/ASE/Aging/cellcounts_scRNA_wgs_genotypes", header=FALSE)
colnames(cell)<- c("ID_SNP", "cell", "chr", "pos", "ref_allele", "alt_allele", "genotype_RNA", "alt", "total", "ID")
cell$cell_id<- paste0(cell$IDnum,"_",cell$cell)

annotation<- read.table("/Users/mharwood/Documents/ASE/scASE/clusterAnnotations", header=FALSE)
colnames(annotation)<- c("cluster_id", "celltype")
annotation$cell_id<- paste0(annotation$cell_id,"-1")
scCountData<- merge(cell, annotation_full, by=c("cell_id"))
scCountData$SNP<- paste0("chr", scCountData$chr, "_", scCountData$pos)

#Filter: 1) total counts > 10 reads and 2) expression data for at least 2 individuals for that SNP
scCountData<- scCountData[scCountData$total>=10,]
scCountData$ID_SNP<- paste0(scCountData$ID,"_",scCountData$SNP)
tt <- table(scCountData$SNP)
scCountData <- subset(scCountData, SNP %in% names(tt[tt>=2]))

#####################
# 	Binomial Test for Adaptive Celltypes	#
#####################

#sum counts if in adaptive cell for the same individual and SNP
adaptive<- scCountData[scCountData$celltype=="CD4+_T"|scCountData$celltype=="CD8+_T"|scCountData$celltype=="B_Cells",]
adaptive<- adaptive %>% group_by(ID_SNP) %>%  dplyr::summarise(across(c(total, der_count), ~sum(.)))
adaptive<- data.frame(adaptive)
colnames(adaptive)<- c("ID_SNP", "total", "der_count")
adaptive$der_count<- as.numeric(adaptive$der_count)
adaptive$total<- as.numeric(adaptive$total)
adaptive<- adaptive[!is.na(adaptive$ID_SNP),]

#binomial test
ngenes <- dim(adaptive)[1]
for (i in 1:ngenes){
  adaptive$pvals[i] <- binom.test(as.numeric(adaptive$der_count[i]),as.numeric(adaptive$total[i]),p=0.5,alternative="two.sided")[["p.value"]]
}

#multiple testing correction
adaptive<-adaptive %>%
  group_by(SNP) %>% 
  mutate(pval.adj = p.adjust (pvals, method='BH'))
sig<- adaptive[adaptive$pval.adj<0.05,]
