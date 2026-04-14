#Run SCT per library and harmony integration across layers

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
#use v5 seurat object and assays.

#updated with adult library.

#read all dataset RDS objects
#dataset paths and names have been partially obsfucated as necessary.
Zhong_dataset <- readRDS("/Zhong_eLife2020_1moMouseBM/Zhong_merged_metaAdded.rds")
Hall_dataset <- readRDS("/Hall_E16E18P0Ad_GSE178951/hall_merged_metaAdded_addAdult_rawSeur.rds")
Ji_dataset <- readRDS("/Ji_CellChemBio2023_8wk20wkMouseBM/Ji_merged_metaAdded_rawSeur.rds")
Song_dataset <- readRDS("/Song_CTRL_GSM6910231/Song_merged_metaAdded_rawSeur.rds")
FanxinLong_dataset <- readRDS("/FanxinLong_10days_18mons/FanxinLong_merged_metaAdded_rawSeur.rds")
Kara_dataset <- readRDS("/Kara_Cell_2023_P4_P14/Kara_P4_P14_rawSeur.rds")

merged_blood_literature <- merge(Zhong_dataset, y = c(Hall_dataset, Ji_dataset, Song_dataset, FanxinLong_dataset, Kara_dataset))
#adds each dataset as a layer, needs to be by orig.ident EACH layer for integration. 

merged_blood_literature <- JoinLayers(merged_blood_literature)
merged_blood_literature[["RNA"]] <- split(merged_blood_literature[["RNA"]], f = merged_blood_literature$orig.ident)


#shared violin plot and then qc filtering, save both full raw data, and filtered rds
VlnPlot(merged_blood_literature, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), raster = FALSE)

saveRDS(merged_blood_literature, file = "rawSeur_merged_blood_literature.rds")


merged_blood_literature_filt <- subset(merged_blood_literature, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nCount_RNA > 400 & nCount_RNA < 60000)

#save
saveRDS(merged_blood_literature_filt, file = "filtSeur_merged_blood_literature.rds")




###Process with SCT and harmony integration


merged_blood_literature_filt <- SCTransform(merged_blood_literature_filt)
merged_blood_literature_filt <- RunPCA(merged_blood_literature_filt)
merged_blood_literature_filt <- RunUMAP(merged_blood_literature_filt, dims = 1:30)
DimPlot(merged_blood_literature_filt, reduction = "umap", group.by = c("orig.ident"))

# integrate datasets
merged_blood_literature_filt <- IntegrateLayers(object = merged_blood_literature_filt, method = HarmonyIntegration, normalization.method = "SCT", verbose = T)
merged_blood_literature_filt <- FindNeighbors(merged_blood_literature_filt, reduction = "harmony", dims = 1:30)
merged_blood_literature_filt <- FindClusters(merged_blood_literature_filt, resolution = 0.8)

merged_blood_literature_filt <- RunUMAP(merged_blood_literature_filt, dims = 1:30, reduction = "harmony")
DimPlot(merged_blood_literature_filt, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

DimPlot(merged_blood_literature_filt, raster = F, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
DimPlot(merged_blood_literature_filt, raster = F, reduction = "umap", group.by = "orig.ident")

#check per library, per mito content, per seq depth clustering. 
FeaturePlot(merged_blood_literature_filt, raster = F, features = "percent.mt")
FeaturePlot(merged_blood_literature_filt, raster = F, features = "nCount_RNA")

saveRDS(merged_blood_literature_filt, file = "blood_wKara_harmony_sctPerLibNorm.rds")
