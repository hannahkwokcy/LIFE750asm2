library(DESeq2)  
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)

#preparing counts data
# read in counts data
counts_data <- read.table('expression/gene_counts.tsv', sep = '\t', header = TRUE)
head(counts_data)
# convert the gene_id column into rownames
rownames(counts_data) <- counts_data$gene_id
counts_data <- counts_data[, -1]

# read in sample info
colData <- read.table('metadata/sample_metadata.tsv', sep = '\t', header = TRUE)
# convert the sample name into rownames
colData <- colData[order(colData$sample), ]
rownames(colData) <- colData$sample

# ensure  row names in colData matches column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# ensure both are in the same order
all(colnames(counts_data) == rownames(colData))

# construct a DESeqDataSet from counts data 
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)
dds

# keep rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# set factor level
dds$condition <- relevel(dds$condition, ref = "normal")

# run deseq
dds <- DESeq(dds)
res <- results(dds)

# order by most significant p-value
res <- res[order(res$padj), ]

# view results
summary(res)

# order by most significant p-value
res0.01 <- results(dds, alpha = 0.01)
res0.01 <- res0.01[order(res0.01$padj), ]
summary(res0.01)

# contrasts
resultsNames(dds)

# MA plot
plotMA(res)

# examine p-values and log2FC distribution in histogram
hist(res$padj, breaks=50, col="grey")
hist(res$log2FoldChange, breaks=50, col="grey")

# rlogTransformatio for heatmap construction
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# set parameters for colors
library(RColorBrewer)
(mycols <-
    brewer.pal(8,"Dark2")[1:length(unique(colData$condition))]
)
# create sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[colData$condition],
          RowSideColors=mycols[colData$condition],
          margin=c(10,10), main="Sample Distance Matrix")

# investigate sample similarities using pca
plotPCA(rld, intgroup=c("condition"))

# find top 5 most upregulated genes
top_upregulated <- res0.01[order(res0.01$log2FoldChange, decreasing = TRUE), ]
top_upregulated <- head(top_upregulated,5)
print(top_upregulated)

# find top 5 most downregulated genes
top_downregulated <- res0.01[order(res0.01$log2FoldChange, decreasing = FALSE), ]
top_downregulated <- head(top_downregulated,5)
print(top_downregulated)

# save list of top 10 genes
top_genes <- c(rownames(top_upregulated), rownames(top_downregulated))

# write results to file
res <- results(dds)
table(res$padj<0.01)
# order by adjusted p-value
res <- res[order(res$padj), ]
# merge with normalized count data
resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds, normalized=TRUE)),
                 by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
# save results to a file
write.table(resdata, file="diffexpr-results.txt",
            sep="\t", quote=F)

# volcano plot to visualise sample distribution
EnhancedVolcano(resdata,
                lab = resdata$Gene,
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-3,3),
                ylim = c(0,32),
                selectLab = top_genes,
                #cutoff values for statistically significant genes
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 1.0,
                max.overlaps = Inf,
                drawConnectors = TRUE,
                labSize = 3.0,
                axisLabSize = 12
)