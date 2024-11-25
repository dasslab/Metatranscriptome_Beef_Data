library(DESeq2)
library(ggplot2)
library(glmpca)
library(plotly)
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library(genefilter)
countData <- read.csv('countdata.csv', header = TRUE, sep = ",", check.names = FALSE)
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Location, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res <- res[order(-res$padj<0.01),]
write.csv(res, "dge.csv")

#PCA Plot
rld <- rlog(dds)
rld <- varianceStabilizingTransformation(dds)
#Plot PCA
pdf("PCA_1.pdf")
pcaData <- plotPCA(rld, intgroup= c("Location", "Plant", "Sample.ID", "name"), ntop=500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Plant, shape=Location, label=name)) +
  geom_point(size=2) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + stat_ellipse(type = "norm", linetype = 2) + theme_bw()
dev.off()



#Heatmap
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
png("heatmap.png",    # create PNG for the heat map
    width = 15*300,        # 5 x 300 pixels
    height = 15*300,
    res = 600,            # 300 pixels per inch
    pointsize = 8) 
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Normal="gray", Tumor="darkgreen")[
             colData(rld)$Location ], margins = c(20,20) )
dev.off()

pdf("MA-Plot.pdf")
plotMA(dds)
dev.off()

pdf("Dispersion.pdf")
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
dev.off()

#Enhanced Volcano
pdf("Enhanced-Volcano.pdf")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5,5),
                ylim = c(0,5),
                pCutoff = 0.05,
                FCcutoff = 2)
dev.off()

#Histrogram
png("Histogram - Padj value.png")
hist( res$padj, breaks=20, col="red" )
dev.off()

png("Histogram - Pvalue.png")
hist( res$pvalue, breaks=20, col="red" )
dev.off()
