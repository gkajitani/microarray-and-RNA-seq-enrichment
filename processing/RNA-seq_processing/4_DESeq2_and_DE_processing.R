############################################################################################
# Using DESeq2 for read count normalization and differential expression analysis - example #
############################################################################################


require("DESeq2")
require(ggplot2)
require (biomaRt)

setwd("~/Desktop/CS_Paper/Todos_datasets_originais/RNASeq")

just.raw.counts = read.csv("read_counts_GSE122736.csv", row.names=1)

head(just.raw.counts)


just.raw.counts <- just.raw.counts[rowSums(just.raw.counts > 10) >= 4, ]


setwd("~/Desktop/CS_Paper/Todos_datasets_originais/RNASeq/Okur_2020")

# Getting only protein coding genes
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

background <- just.raw.counts

background$Genes <- rownames(background)

background <- merge(genes, background, by.x = "external_gene_name", by.y = "Genes")

write.csv(background, "background_Okur_2020.csv")


setwd("~/Desktop/CS_Paper/Todos_datasets_originais/RNASeq/Okur_2020")


meta.data = read.csv(file="meta.data_Okur_2020.csv", row.names = 1)

head(meta.data)

count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ condition)

count.data.set.object <- DESeq(count.data.set)

## estimating size factors

vsd <- vst(count.data.set.object)
norm.data = assay(vsd)

head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt", row.names=TRUE,col.names=NA,quote=FALSE)
sampleDists <- dist(t(norm.data),  method = "euclidean")

reversed_rows_columns = (t(norm.data))
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("condition")) 


# Comparison (note: reference data is always the first one in the "contrast" option - i.e. positive fold change indicates upregulation of gene expression in the reference)
res <- results(count.data.set.object, contrast=c("condition", "CSB", "WT"))

csb <- as.data.frame(res)

write.csv(csb,"all_CS1AN.csv")

res
summary(res)
res = na.omit(res)
head(res)
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
summary(resfiltered)
sum(resfiltered$padj > 0.05)

res.filtered.ordered = resfiltered[order(resfiltered$padj),] 
head(res.filtered.ordered)
write.table(res.filtered.ordered, sep="\t",file="Results_filtered_ordered.txt", row.names=TRUE,col.names=NA,quote=FALSE)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)

resOrdFilt.data.frame$RowNames <- rownames(resOrdFilt.data.frame)


merged_data <- merge(genes, resOrdFilt.data.frame, by.x = "external_gene_name", by.y = "RowNames")

write.csv(merged_data,"Okur_todos.csv")


resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)
head(resOrdFilt.data.frame)
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram()
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-10,10))
resOrdFiltLFC = subset(resOrdFilt.data.frame,log2FoldChange >0.5 | log2FoldChange < -0.5)
dim(resOrdFiltLFC)
ggplot(resOrdFiltLFC, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-10,10))
write.table(resOrdFiltLFC, sep="\t",file="Results_filtered_ordered_LFC.txt", row.names=TRUE,col.names=NA,quote=FALSE)
summary(res.filtered.ordered)

resOrdFiltLFC$RowNames <- rownames(resOrdFiltLFC)


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)


merged_data <- merge(genes, resOrdFiltLFC, by.x = "external_gene_name", by.y = "RowNames")


UpInCSB = subset(merged_data, log2FoldChange >0.5)
DownInCSB = subset(merged_data, log2FoldChange <(-0.5))



write.csv(merged_data, "Okur_2020_LFC.csv")


write.csv(UpInCSB, "UpInCSB_Okur.csv")
write.csv(DownInCSB, "DownInCSB_Okur.csv")

