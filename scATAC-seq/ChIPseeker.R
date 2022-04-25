library(ChIPseeker)
library(org.At.tair.db)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(clusterProfiler)
library(ggupset)

txdb <- TxDb.Athaliana.BioMart.plantsmart28
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
gis_1<-readPeakFile("1.peaks.narrowPeak")


gis_1@seqnames@values <- promoter@seqnames@values
gis_1@seqinfo@seqnames <- promoter@seqinfo@seqnames

tagMatrix_1 <- getTagMatrix(gis_1, windows=promoter)
#plotAvgProf(tagMatrix, xlim=c(-1000, 1000),
 #           conf=0.95,resample = 1000,
  #          xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")


peakAnno_1 <- annotatePeak(gis_1, tssRegion=c(-1000, 1000),TxDb=txdb, 
                         annoDb="org.At.tair.db",level="gene")
plotAnnoPie(peakAnno_1)


upsetplot(peakAnno_1)
df <- as.data.frame(peakAnno_1)
gene_1 <- df[,19]
gene_1 <- unique(gene_1)
up_geo_1 <- clusterProfiler::enrichGO(gene = gene_1, keyType = "TAIR",
                                    OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                    ont = "BP",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)
enrichplot::dotplot(up_geo_1, showCategory =30)


gis_6<-readPeakFile("6.peaks.narrowPeak")

gis_6@seqnames@values <- promoter@seqnames@values
gis_6@seqinfo@seqnames <- promoter@seqinfo@seqnames
tagMatrix_6 <- getTagMatrix(gis_6, windows=promoter)
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
 #           conf=0.95,resample = 1000,
  #          xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

peakAnno_6 <- annotatePeak(gis_6, tssRegion=c(-1000, 1000),TxDb=txdb, 
                         annoDb="org.At.tair.db",level="gene")


plotAnnoPie(peakAnno_6)
upsetplot(peakAnno_6)
df <- as.data.frame(peakAnno_6)
gene_6 <- df[,19]
gene_6 <- unique(gene_6)
up_geo <- clusterProfiler::enrichGO(gene = gene_6, keyType = "TAIR",
                                    OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                    ont = "BP",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)
enrichplot::dotplot(up_geo_6, showCategory =30)

up_geo <- clusterProfiler::enrichGO(gene = only_1, keyType = "TAIR", 
                                    OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                    ont = "BP",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(up_geo, showCategory =10)

up_geo_6 <- clusterProfiler::enrichGO(gene = only_6, keyType = "TAIR", 
                                    OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                    ont = "All",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(up_geo_6, showCategory =10)
