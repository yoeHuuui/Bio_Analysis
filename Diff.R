library(ggpubr)

DESeq <- read.table("DESeq2_logFC.csv", sep = ',')
DESeq <- DESeq[-1, ]
DESeq <- DESeq[, c(1, 3, 7)]
colnames(DESeq) <- c("Gene", "log2FC", "pvalue")

subset(DESeq,pvalue < 0.01) -> DESeq_diff
subset(DESeq_diff,log2FC < -2) -> DESeq_down
subset(DESeq_diff,log2FC > 2) -> DESeq_up

Gfold <- read.table("data.diff")
colnames(Gfold) <- c("Symbol", "Gene", "Gvalue", "E-FDR", "log2FC", "1", "2")
Gfold <- Gfold[, c(2, 3, 4, 5)]

subset(Gfold,abs(Gvalue) > 0) -> Gfold_diff
subset(Gfold_diff,log2FC < -2) -> Gfold_down
subset(Gfold_diff,log2FC > 2) -> Gfold_up

common_diff <- intersect(x=DESeq_diff$Gene, y=Gfold_diff$Gene)

Gfold_fc <- Gfold_diff[Gfold_diff$Gene %in% common_diff,]$log2FC
DESeq_fc <- DESeq_diff$log2FC

FC_data <- data.frame(-Gfold_fc, as.numeric(DESeq_fc))
colnames(FC_data) <- c("G", "D")

ggscatter(data=FC_data, x="G", y="D", add = 'reg.line', conf.int = TRUE,
          add.params = list(color = "#85483650", size = 2), cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = -5), 
          color = "#2F568875", size = 3)+
  xlab("Gfold_FC")+ylab("DESeq_FC")
ggsave("Pearson.png", dpi = 600, width = 8, height = 6)  
