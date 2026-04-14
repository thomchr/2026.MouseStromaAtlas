##Extract the SCVI embeddings and add to the fully processed harmony-integrated seurat object
#downstream normalization and visualization performed in R computing framework and Seurat packages


library(zellkonverter)
library(Seurat)

sce <- readH5AD("./SCVI_integration/scvi_integrated_cellranger_v4.h5ad")

# Convert to Seurat
scvi_integrated_seur <- as.Seurat(sce, counts = NULL, data = 'X')

#confirm concordance of metadata and embeddings
#confirm concordance of dimplot

#saveRDS(seurat_obj, file = 'scvi_integrated_cellranger_epoch400_v4_seur.rds')


#readin both seurat objects
#scvi_integrated_seur <- readRDS("scvi_integrated_cellranger_epoch400_v4_seur.rds")

#harmony integrated seurat
harmony_integrated_seur <- readRDS("blood_wKara_harmony_sctPerLibNorm.rds")
#confirm default assay is 'RNA' or set with DefaultAssay()

#add reductions, leiden clusters, Xscvi_batch and xscvi_labels
#might need to subset both objects for only shared
#modify the rownames of the scvi-integrated meta to match the harmony seurat naming
write.csv(harmony_integrated_seur[[]], file = 'harmony_meta_naming.csv')
write.csv(scvi_integrated_seur[[]], file = 'scvi_integ_meta_naming.csv')

#change names in the scvi object, all 18+1 libraries

new_names <- gsub("-1-Zhong_3M", "-1_3_1", colnames(scvi_integrated_seur))
new_names <- gsub("-1-Zhong_1M", "-1_1_1", new_names)
new_names <- gsub("-1-Zhong_1_5M", "-1_2_1", new_names)
new_names <- gsub("-1-Zhong_16M", "-1_4_1", new_names)


new_names <- gsub("-1-Ji_WT_8wks", "-1_1_3", new_names)
new_names <- gsub("-1-Ji_WT_20wks", "-1_2_3", new_names)


new_names <- gsub("-1-CTRL_BMStroma_Song", "-1_4", new_names)

new_names <- gsub("-1-FL_neonatal_10_days", "-1_2_5", new_names)
new_names <- gsub("-1-FL_aging_18_mons", "-1_1_5", new_names)


new_names <- gsub("-1-P0_CD45negTer119neg_BMcells_comb", "-1_6_2", new_names)
new_names <- gsub("-1-P0_lineageKitpos_BMcells", "-1_7_2", new_names)
new_names <- gsub("-1-emb_day16_5_lineageKitpos_BMcells2", "-1_3_2", new_names)
new_names <- gsub("-1-emb_day18_5_lineageKitpos_BMcells", "-1_5_2", new_names)
new_names <- gsub("-1-emb_day18_5_CD45negTer119neg_BMcells", "-1_4_2", new_names)
new_names <- gsub("-1-emb_day16_5_CD45negTer119neg_BMcells", "-1_1_2", new_names)
new_names <- gsub("-1-emb_day16_5_lineageKitpos_BMcells1", "-1_2_2", new_names)
new_names <- gsub("-1-Adult_LineageNegKitPos_BMCells", "-1_8_2", new_names)


new_names <- gsub("-1-Kara_P4", "-1_1_6", new_names)
new_names <- gsub("-1-Kara_P14", "-1_2_6", new_names)

#write.csv(as.data.frame(new_names), file = "new_names_check_meta.csv")

colnames(scvi_integrated_seur) <- new_names

#now can add into harmony object
#maybe add BOTH X_umap and X_scVI to new harmony object
Reductions(scvi_integrated_seur)
emb_scvi <- Embeddings(scvi_integrated_seur[["X_scVI"]])
head(emb_scvi)
# Get common cells
common_cells <- intersect(rownames(harmony_integrated_seur@meta.data), rownames(emb_scvi))
#169,020

# Subset embedding to only those cells
emb_aligned <- emb_scvi[common_cells, ]
harmony_integrated_seur[["scvi_embed"]] <- CreateDimReducObject(
  embeddings = emb_aligned,
  key = "scvi_",
  assay = DefaultAssay(harmony_integrated_seur)
)

#add this and meta back to back harmony seurat and save
harmony_integrated_seur[['scvi_leiden_1']] <- scvi_integrated_seur$leiden

#add the x_umap
emb_umap_scvi <- Embeddings(scvi_integrated_seur[["X_umap"]])
head(emb_umap_scvi)
# Get common cells
common_cells <- intersect(rownames(harmony_integrated_seur@meta.data), rownames(emb_umap_scvi))
#169,020

# Subset embedding to only those cells
emb_aligned <- emb_umap_scvi[common_cells, ]
harmony_integrated_seur[["scvi_umap"]] <- CreateDimReducObject(
  embeddings = emb_aligned,
  key = "scviUMAP_",
  assay = DefaultAssay(harmony_integrated_seur)
)

DimPlot(harmony_integrated_seur, raster = F, reduction = 'scvi_embed', group.by = 'scvi_leiden_1')
p1 <- DimPlot(harmony_integrated_seur, raster = F, reduction = 'scvi_umap', group.by = 'scvi_leiden_1', label = T) + NoLegend()
p2 <- DimPlot(harmony_integrated_seur, raster = F, reduction = 'umap', group.by = 'seurat_clusters', label = T) + NoLegend()
p1 + p2

p1 <- DimPlot(harmony_integrated_seur, raster = F, reduction = 'scvi_umap', group.by = 'orig.ident', label = F) + NoLegend()
p2 <- DimPlot(harmony_integrated_seur, raster = F, reduction = 'umap', group.by = 'orig.ident', label = F) + NoLegend()
p1 + p2

saveRDS(harmony_integrated_seur, file = 'harmony_integrated_scvi.rds')