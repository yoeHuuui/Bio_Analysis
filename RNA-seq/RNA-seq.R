library(DESeq2)
library(ggplot2)
library(apeglm)
library(gplots)
library(ggpubr)
library(ggthemes)
library("org.At.tair.db", character.only = TRUE)

Matrix1 <- read.table("matrix1.out", header=TRUE, row.names=1)
Matrix2 <- read.table("matrix2.out", header=TRUE, row.names=1)
colnames(Matrix1) <- c("Whole")
colnames(Matrix2) <- c("Root")
Data <- cbind(Matrix2, Matrix2, Matrix1, Matrix1)
colnames(Data) <- c("Root_1", "Root_2", "Whole_1", "Whole_2") 
condition <- factor(c("Root", "Root", "Whole", "Whole"))
coldata <- data.frame(row.names=colnames(Data), condition)

dds <- DESeq2::DESeqDataSetFromMatrix(countData=Data, 
                              colData=coldata, design= ~condition)

dds <- dds[rowSums(BiocGenerics::counts(dds)) > 100, ] 
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds) 
resdata <- data.frame(results(dds, lfcThreshold=1, alpha=0.05)) 
write.csv(resdata,file= "DESeq2_logFC.csv")


threshold <- as.factor(ifelse(resdata$padj < 0.001 & 
                                
                                abs(resdata$log2FoldChange) >= 2 ,
                              
                              ifelse(resdata$log2FoldChange >= 2 ,
                                     'Up','Down'),'Not'))

ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),colour=threshold)) +
  xlab("log2(Fold Change)")+ylab("-log10(qvalue)") +
  geom_point() +
  ylim(0,50) + xlim(-12,12) +
  scale_color_manual(values=c("#2F5688", "#BBBBBB", "#CC0000"))+theme_base()+
  geom_hline(yintercept = -log10(0.001), linetype = "dashed")+
  geom_vline(xintercept = c(-2, 2), linetype = "dashed")
ggsave("Vocalno_DE.png", dpi = 600, width = 8, height = 6)


subset(resdata,pvalue < 0.001) -> diff
subset(diff,log2FoldChange < -2) -> down
subset(diff,log2FoldChange > 2) -> up

down_geo <- clusterProfiler::enrichGO(gene = rownames(down), keyType = "TAIR",
                                    OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                    ont = "BP",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(down_geo, showCategory =30)
ggsave("Down_GO.png", dpi = 600, width = 8, height = 6)

up_geo <- clusterProfiler::enrichGO(gene = rownames(up), keyType = "TAIR",
                                      OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                      ont = "BP",
                                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(up_geo, showCategory =30)
ggsave("Up_GO.png", dpi = 600, width = 8, height = 6)

