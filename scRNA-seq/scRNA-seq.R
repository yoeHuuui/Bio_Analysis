library(dplyr)
library(Seurat)
library(monocle3)
library(ggplot2)
library("org.At.tair.db", character.only = TRUE)


Ath.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
Ath <- CreateSeuratObject(counts = Ath.data, project = "Ath", 
                           min.cells = 10, min.features = 300)
Ath[["nCercent_mt"]] <- PercentageFeatureSet(Ath, pattern = "*MG")
VlnPlot(Ath, features = c("nFeature_RNA", "nCount_RNA", "nCercent_mt"), 
        ncol = 3)
ggsave("QC.svg", width = 8, 
       height = 6)

Ath <- subset(Ath, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & nCercent_mt < 5)
Ath <- NormalizeData(Ath, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

Ath <- FindVariableFeatures(Ath, selection.method = "vst", 
                            nfeatures = 2000)
VariableFeaturePlot(Ath)
ggsave("VariableFeature.svg", width = 8, 
       height = 6, dpi = 1200)

all.genes <- rownames(Ath)
Ath <- ScaleData(Ath, features = all.genes)
Ath <- RunPCA(Ath, features = VariableFeatures(object = Ath), 
              npcs = 30)
DimPlot(Ath, reduction = "pca")
ggsave("Dim.svg", width = 8, 
       height = 6)

#DimHeatmap(Ath, dims = 1:6, cells = 1000, balanced = TRUE)

Ath <- JackStraw(Ath, num.replicate = 10)
Ath <- ScoreJackStraw(Ath, dims = 1:20)
JackStrawPlot(Ath, dims = 1:20)
ElbowPlot(Ath)
ggsave("Elbow.svg", width = 8, 
       height = 6)

Ath <- FindNeighbors(Ath, dims = 1:20)
Ath <- FindClusters(Ath, resolution = 0.7)

Ath <- RunUMAP(Ath, dims = 1:10)
Ath <- RunTSNE(Ath, dims = 1:10)

DimPlot(Ath, reduction = "umap")
ggsave("Umap.svg", width = 8, 
       height = 6)
DimPlot(Ath, reduction = "tsne")
ggsave("Tsne.svg", width = 8, 
       height = 6)

Ath.markers <- FindAllMarkers(Ath, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)

top10 <- Ath.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Ath, features = top10$gene) + NoLegend()
ggsave("DoHeatmap.png", width = 18, height = 24, dpi=300)

features <- c("AT2G16750", "AT4G30170")
RidgePlot(Ath, features = features, ncol = 2)
ggsave("Ridge.svg", width = 8, 
       height = 6)

saveRDS(Ath, file = "Ath.rds")


data <- as(as.matrix(Ath@assays$RNA@counts), 'sparseMatrix')
pd <-  Ath@meta.data
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)
cds = preprocess_cds(cds, num_dim = 30)
plot_pc_variance_explained(cds)
cds = reduce_dimension(cds, reduction_method="UMAP")
plot_cells(cds, color_cells_by="seurat_clusters")
ggsave("Umap_mo.svg", width = 8, 
       height = 6)

cds = reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")
ggsave("Tsne_mo.svg", width = 8, 
       height = 6)

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds)


plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=2)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=partitions(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)



Top_gene <- clusterProfiler::enrichGO(top10$gene, keyType = "TAIR",
                                      OrgDb = BiocGenerics::get("org.At.tair.db"), 
                                      ont = "BP",
                                      pAdjustMethod = "BH", qvalueCutoff = 0.05)
enrichplot::dotplot(Top_gene, showCategory =30)


