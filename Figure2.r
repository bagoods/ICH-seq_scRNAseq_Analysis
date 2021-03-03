library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library('data.table')
library(stringr)
library(dplyr)
library(RColorBrewer)


library("ComplexHeatmap")
library("ggplot2")

library("ggbiplot")

project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
project
workspace
bucket

list.files(path = ".")

# use below path for first run through analysis, the second for subsequent ones
#RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/*.rds"

RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/*.rds"
system(sprintf("gsutil ls %s", RDS_Files), intern=T)

system("mkdir $(pwd)/RDS_Files")
system(sprintf("gsutil -m cp %s $(pwd)/RDS_Files/", RDS_Files))

list.files(path = "./RDS_Files")

#seurat.all.Tcells <- readRDS("./RDS_Files/seuratobj_alldata_filtered_Tcells.rds") #import for re-analysis originally 

seurat.all.data <- readRDS("./RDS_Files/seuratobj_seurat_all_data_TcellFigure_correct.rds") #import post module scoring 


#seurat.all.Tcells

# subset to remove NK cells (only on first pass)
#Idents(object=seurat.all.Tcells) <- "idents_res2_consensus"
#seurat.all.data <- subset(seurat.all.Tcells, idents = c("activated T cells", "proliferating T cells", "CD8 T cells", "heat shock response T cells", 
#"T cells 4", "T cells 2", "ribosomal-enriched T cells", "XIST+ T cells", 
#"T cells 1", "T cells 5", "T cells 6", "T cells 7"))
#saveRDS(seurat.all.data,file="seuratobj_seurat_all_data_TcellFigure_correct.rds")
#system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure_correct.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



seurat.all.data

unique(seurat.all.data@meta.data$idents_res2_consensus)

# colors 
cols.hours = c("lightsteelblue1","darkslategray1","cadetblue2","dodgerblue2","blue1","darkblue","gold")
cols.compartment = c("brown","darkgray")
cols.platform = c("rosybrown","peru")

## normalize and initial clustering 
seurat.all.data <- NormalizeData(seurat.all.data, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.all.data <- FindVariableFeatures(seurat.all.data, selection.method = "vst", nfeatures = 2000)

# Identify the 30 most highly variable genes
top_merge <- head(VariableFeatures(seurat.all.data), 30)

# plot variable features with and without labels
plotmerge1 <- VariableFeaturePlot(seurat.all.data)
plotmerge2 <- LabelPoints(plot = plotmerge1, points = top_merge, repel = TRUE)
plotmerge2
ggsave("Variable_Genes_Merge.pdf", useDingbats = F)
system(paste0("gsutil cp Variable_Genes_Merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# scale and PCA, can regress stuff out here too (do not regress for now)
seurat.all.data <- ScaleData(seurat.all.data, features = rownames(seurat.all.data))
seurat.all.data <- RunPCA(seurat.all.data, features = VariableFeatures(object = seurat.all.data))

ElbowPlot(seurat.all.data)
ggsave("Elbow_merge.pdf", useDingbats = F)
system(paste0("gsutil cp Elbow_merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


DimHeatmap(seurat.all.data, dims = 1:20, cells = 500, balanced = TRUE)
ggsave("PC_heatmaps.pdf", useDingbats = F, height = 10, width = 10, units = "in")
system(paste0("gsutil cp PC_heatmaps.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


head(seurat.all.data@meta.data)

Idents(object=seurat.all.data) <- "ArrayName"

DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_array.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_array.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "Platform"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_platform.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_platform.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# PCA plots colored by compartment
Idents(object=seurat.all.data) <- "Compartment"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# PCA plots colored by SNN2 cluster 
Idents(object=seurat.all.data) <- "idents_res2_consensus"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_consensus.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_consensus.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_consensus.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_consensus.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# PCA plots colored by time
Idents(object=seurat.all.data) <- "Hours"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_hours.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_hours.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


## Jackstraw sig PCs (takes a really long time to run - all 20 PCs are sig!)
seurat.all.data <- JackStraw(seurat.all.data, num.replicate = 100)
seurat.all.data <- ScoreJackStraw(seurat.all.data, dims = 1:20)
JackStrawPlot(seurat.all.data, dims = 1:20)
ggsave("JackstrawPlot.pdf", useDingbats = F)
system(paste0("gsutil cp JackstrawPlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_seurat_all_data_TcellFigure_correct.rds")
system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure_correct.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# elbow 
ElbowPlot(seurat.all.data)


seurat.all.data <- RunTSNE(seurat.all.data, reduction = "pca", dims = 1:18)
seurat.all.data <- FindNeighbors(seurat.all.data, reduction = "pca", dims = 1:18)


# QC feature plots
FeaturePlot(seurat.all.data, features = c("percent.mt"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_percentmito.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_percentmito.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# colored by other features
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform)
ggsave("UMAP_platform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform")
ggsave("UMAP_platform_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment)
ggsave("UMAP_compartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment")
ggsave("UMAP_compartment_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8)
ggsave("UMAP_arrayname.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours)
ggsave("UMAP_hours.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours")
ggsave("UMAP_hours_splitbyhrs.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate")
ggsave("UMAP_hours_splitbyreplicate.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform")
ggsave("UMAP_hours_splitbyplatform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment")
ggsave("UMAP_hours_splitbycompartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# cell ID
DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8)
ggsave("UMAP_cellID.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Hours")
ggsave("UMAP_cellID_splitbyhrs.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyhrs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Replicate")
ggsave("UMAP_cellID_splitbyreplicate.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Platform")
ggsave("UMAP_cellID_splitbyplatform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Compartment")
ggsave("UMAP_cellID_splitbycompartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbycompartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



# colored by other features (no legend)
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform) + NoLegend()
ggsave("UMAP_platform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform") + NoLegend()
ggsave("UMAP_platform_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment) + NoLegend()
ggsave("UMAP_compartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment") + NoLegend()
ggsave("UMAP_compartment_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8) + NoLegend()
ggsave("UMAP_arrayname_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours) + NoLegend()
ggsave("UMAP_hours_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours") + NoLegend()
ggsave("UMAP_hours_splitbyhrs_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate") + NoLegend()
ggsave("UMAP_hours_splitbyreplicate_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform") + NoLegend()
ggsave("UMAP_hours_splitbyplatform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment") + NoLegend()
ggsave("UMAP_hours_splitbycompartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# cell ID
DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8)+ NoLegend()
ggsave("UMAP_cellID_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Hours")+ NoLegend()
ggsave("UMAP_cellID_splitbyhrs_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyhrs_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Replicate")+ NoLegend()
ggsave("UMAP_cellID_splitbyreplicate_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Platform")+ NoLegend()
ggsave("UMAP_cellID_splitbyplatform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Compartment")+ NoLegend()
ggsave("UMAP_cellID_splitbycompartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbycompartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



Idents(object=seurat.all.data) <- "ArrayName"
DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5)

cols.array.pt.ctl = c("lightsteelblue1","lightsteelblue1",
                     "darkslategray1","darkslategray1",
                     "cadetblue2","cadetblue2","cadetblue2",
                     "dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2",
                     "blue1","blue1","blue1","blue1",
                     "darkblue","darkblue","darkblue",
                     "gold","gold","gold","gold",
                     "orange")

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.array.pt.ctl)


# tSNEs colored with the PT and controls 
DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.array.pt.ctl)  
ggsave("tsne_time_collapse_PThighlight.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_time_collapse_PThighlight.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.array.pt.ctl) + NoLegend()
ggsave("tsne_time_collapse_PThighlight_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_time_collapse_PThighlight_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)





# cluster IDs (0.4, 0.5, 0.6, 0.8) #0.6 is the one 
seurat.all.data <- FindClusters(seurat.all.data, resolution = 0.6)


colnames(seurat.all.data@meta.data)

unique(seurat.all.data@meta.data$RNA_snn_res.0.6)

#10
#colors.cluster.new <- c("indianred1","hotpink2","darkgoldenrod1","olivedrab3","springgreen4","lightseagreen",
#                        "deepskyblue","royalblue","deepskyblue","darkorchid4")

#12
colors.cluster.new <- c("antiquewhite3","indianred1","hotpink2","darkgoldenrod","springgreen4","darkolivegreen",
                        "lightseagreen","deepskyblue","royalblue","blue4","purple","darkorchid4") 

#14
#colors.cluster.new <- c("antiquewhite3","indianred1","hotpink2","gold3","darkgoldenrod","darkseagreen","springgreen4",
#                        "darkolivegreen","lightseagreen","deepskyblue","royalblue","blue4","purple","darkorchid4")
    

# all together
DimPlot(seurat.all.data, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colors.cluster.new)  
ggsave("tsne_clusterID_res06.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colors.cluster.new) + NoLegend()
ggsave("tsne_clusterID_res06_no_legend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06_no_legend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#split by compartment
DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment", cols = colors.cluster.new)  
ggsave("tsne_clusterID_res06_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment", cols = colors.cluster.new) + NoLegend()
ggsave("tsne_clusterID_res06_no_legend_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06_no_legend_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



colnames(seurat.all.data@meta.data)

## Try 0.4, 0.6, 0.8
Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
seurat.all.data.markers <- FindAllMarkers(seurat.all.data, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5)

write.table(seurat.all.data.markers,file="seurat.all.markers.wilcox.res06.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.markers.wilcox.res06.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


head(seurat.all.data.markers)

#seurat.all.data.markers %>% group_by(cluster) %>% top_n(10, p_val_adj) -> top10

top10 <- seurat.all.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


top10

DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_res06.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

#DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))
#ggsave("Exp_Heatmap_res04_BluYR.pdf", useDingbats = F)
#system(paste0("gsutil cp Exp_Heatmap_res04_BluYR.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


seurat.all.data.averages <- AverageExpression(seurat.all.data, return.seurat = TRUE)

DoHeatmap(seurat.all.data.averages,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_Avg_res06.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_Avg_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


seurat.all.data <- BuildClusterTree(seurat.all.data, verbose = FALSE, reorder = TRUE)
tree <- PlotClusterTree(seurat.all.data)
ggsave("ClusterTree_res06.pdf", useDingbats = F)
system(paste0("gsutil cp ClusterTree_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_seurat_all_data_TcellFigure_correct.rds")
system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure_correct.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("TRAC", "CD4", "CD8A", "IL7R",
                                          "TRBC2","CD3E","CD3G",
                                         "GNLY","NKG7","GZMA"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_Tcell_FeaturePlot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Tcell_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



FeaturePlot(seurat.all.data, features = c("PTPRC","CD8A","CD4","CCR7"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("Feature_plot_reveiwercomments.pdf", useDingbats = F)
system(paste0("gsutil cp Feature_plot_reveiwercomments.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# dotplot of receptor and ligands 

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
DotPlot(seurat.all.data, features = c("CCR9","CCR7","CXCR5","CXCR4","CCR10","CX3CR1","CCR8",
                                     "CCR5","CCR3","CCR2","CCR1","CXCR6","CXCR3","CCL1","CCL2","CCL18","MAP2K1",
                                     "CCL4","CCL7","CCL5","CXCL16","CXCL10")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("DotPlot.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp DotPlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("CD4", "CD8A", "CCR7","SELL", "PTPRC", "IL7R"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("Feature_plot_reveiwercomments_2.pdf", useDingbats = F)
system(paste0("gsutil cp Feature_plot_reveiwercomments_2.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"

VlnPlot(seurat.all.data, features = c(features = c("CD4", "CD8A", "CCR7","SELL", "PTPRC", "IL7R")), pt.size = 0.1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_Reviewer3_res06.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_Reviewer3_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c(features = c("TRAC", "CD4", "CD8A", "IL7R",
                                          "TRBC2","CD3E","CD3G",
                                         "GNLY","NKG7","GZMA")), pt.size = 0, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c(features = c("TRAC", "CD4", "CD8A", "IL7R",
                                          "TRBC2","CD3E","CD3G",
                                         "GNLY","NKG7","GZMA")), pt.size = 0.1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_points_res06.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_points_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# violin plots of marker genes
VlnPlot(seurat.all.data, features = c("TXNIP","PABPC1","S100A8","BHLHE40","GZMA",
"HSPA1A","TIMP1","MALAT1","ZNF90","HLA-DQA1","STMN1","RP11-297H3.3"), pt.size = 0, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_markers.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_markers.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)




unique(colnames(seurat.all.data@meta.data))

## stacked barplot of clusters by donor, compartment, platform, and hours 
# donor (replicate)
ggplot(seurat.all.data@meta.data, aes(x=RNA_snn_res.0.6, fill=Replicate)) + geom_bar(position = "fill")
ggsave("StackedBarClusters_res06_replicate.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# compartment
ggplot(seurat.all.data@meta.data, aes(x=RNA_snn_res.0.6, fill=Compartment)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.compartment)
ggsave("StackedBarClusters_res06_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# platform
ggplot(seurat.all.data@meta.data, aes(x=RNA_snn_res.0.6, fill=Platform)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.platform)
ggsave("StackedBarClusters_res06_platform.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# hours
ggplot(seurat.all.data@meta.data, aes(x=RNA_snn_res.0.6, fill=Hours)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.hours)
ggsave("StackedBarClusters_res06_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

# proper colors over array colored by consensus ID
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("StackedBarClusters_CellID_Array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_CellID_Array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)




# proper colors over array colored by cluster number
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "Compartment"
suerat.all.data.subset.BLD <- subset(seurat.all.data, idents = c("BLD"))
suerat.all.data.subset.HEF <- subset(seurat.all.data, idents = c("HEF"))


colors.cluster.new.BLD <- c("indianred1","hotpink2","darkgoldenrod1","olivedrab3","springgreen4","lightseagreen",
                        "royalblue","darkorchid4")

ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Array_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Hours_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Hours_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Array_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


unique(seurat.all.data@meta.data$ArrayName)

Idents(object=seurat.all.data) <- "ArrayName"
suerat.all.data.subset.BLD.CT <- subset(seurat.all.data, idents = c("DayF-500-BLD-PT-Nova"))
suerat.all.data.subset.BLD.PT <- subset(seurat.all.data, idents = c("DayF-500-BLD-00-Nova","DayF-500-BLD-01-Nova",
                                                                   "DayF-500-BLD-02-Nova","DayF-500-BLD-03-Nova"))



colors.cluster.new.BLD.CT <- c("antiquewhite3","hotpink2","darkgoldenrod","springgreen4","blue4","purple") 

ggplot(suerat.all.data.subset.BLD.CT@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.BLD.CT@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.BLD.CT)
ggsave("StackedBarClusters_ClusterID_Array_BLD_CT.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_BLD_CT.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#12
colors.cluster.new <- c("antiquewhite3","indianred1","hotpink2","darkgoldenrod","springgreen4","darkolivegreen",
                        "lightseagreen","deepskyblue","royalblue","blue4","purple","darkorchid4") 

colors.cluster.new.BLD.PT <- c("antiquewhite3","indianred1","hotpink2","darkgoldenrod",
                               "springgreen4","deepskyblue","blue4","purple")


ggplot(suerat.all.data.subset.BLD.PT@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.BLD.PT@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.BLD.PT)
ggsave("StackedBarClusters_ClusterID_Array_BLD_PT.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_BLD_PT.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


library(VISION)

#load .gmt file from bucket 
GMT_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/*.gmt"
system(sprintf("gsutil ls %s", GMT_Files), intern=T)

system("mkdir $(pwd)/GMT_Files")
system(sprintf("gsutil -m cp %s $(pwd)/GMT_Files/", RDS_Files))

list.files(path = "./GMT_Files")

signatures <- c("GMT_Files/c7.all.v7.1.symbols.gmt")

vision.obj <- Vision(seurat.all.data, signatures = signatures)


vision.obj <- analyze(vision.obj)

viewResults(vision.obj)

colnames(seurat.all.data@meta.data)

seurat.all.data@meta.data$Cano_Gamez_Orange1 <- NULL
seurat.all.data@meta.data$Cano_Gamez_Green1 <- NULL
seurat.all.data@meta.data$Cano_Gamez_Green_High1 <- NULL
seurat.all.data@meta.data$Cano_Gamez_Green_Low1 <- NULL


gene.sets <- read.csv(pipe('gsutil cat gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/TcellSignatures.csv'))


head(gene.sets)
unique(colnames(gene.sets))

## add in additional cano_gamez scores 

Cano_Gamez_Orange_Resting <- list(gene.sets$Cano_Gamez_Orange_Resting)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Cano_Gamez_Orange_Resting, name = "Cano_Gamez_Orange_Resting", seed.use = 1988)

Cano_Gamez_Green_Resting <- list(gene.sets$Cano_Gamez_Green_Resting)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Cano_Gamez_Green_Resting, name = "Cano_Gamez_Green_Resting", seed.use = 1988)

Cano_Gamez_Orange <- list(gene.sets$Cano_Gamez_Orange)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Cano_Gamez_Orange, name = "Cano_Gamez_Orange", seed.use = 1988)

Cano_Gamez_Blue <- list(gene.sets$Cano_Gamez_Blue)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Cano_Gamez_Blue, name = "Cano_Gamez_Blue", seed.use = 1988)

Cano_Gamez_Green <- list(gene.sets$Cano_Gamez_Green)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Cano_Gamez_Green, name = "Cano_Gamez_Green", seed.use = 1988)


Szabo_Treg <- list(gene.sets$Szabo_Treg)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_Treg, name = "Szabo_Treg", seed.use = 1988)

Szabo_Treg <- list(gene.sets$Szabo_Treg)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_Treg, name = "Szabo_Treg", seed.use = 1988)

Szabo_CD4NV.CM <- list(gene.sets$Szabo_CD4NV.CM)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_CD4NV.CM, name = "Szabo_CD4NV.CM", seed.use = 1988)

Szabo_CD4.CD8 <- list(gene.sets$Szabo_CD4.CD8)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_CD4.CD8, name = "Szabo_CD4.CD8", seed.use = 1988)

Szabo_IFNResponse <- list(gene.sets$Szabo_IFNResponse)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_IFNResponse, name = "Szabo_IFNResponse", seed.use = 1988)

Szabo_Proliferation <- list(gene.sets$Szabo_Proliferation)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_Proliferation, name = "Szabo_Proliferation", seed.use = 1988)

Szabo_CD8Cytotoxic <- list(gene.sets$Szabo_CD8Cytotoxic)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_CD8Cytotoxic, name = "Szabo_CD8Cytotoxic", seed.use = 1988)

Szabo_CD8Cytokine <- list(gene.sets$Szabo_CD8Cytokine)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Szabo_CD8Cytokine, name = "Szabo_CD8Cytokine", seed.use = 1988)


colnames(seurat.all.data@meta.data)

#produce.module.score = function(gene.sets, seurat.object){
# new.seurat.object <- seurat.object
# for (i in colnames(gene.sets)){
#   list.genes <- list(unique(gene.sets[[i]]))
#   new.seurat.object <- AddModuleScore(object = new.seurat.object, features = list.genes, name = i, seed = 1996)
# }
# return(new.seurat.object)
# saveRDS(new.seurat.object,file="seuratobj_seurat_all_data_TcellFigure_modules.rds")
# system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure_modules.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)
#}

FeaturePlot(seurat.all.data, features = c("Szabo_Treg1","Szabo_CD4NV.CM1","Szabo_CD4.CD81","Szabo_IFNResponse1","Szabo_Proliferation1",
                                         "Szabo_CD8Cytotoxic1","Szabo_CD8Cytokine1"), pt.size = 0.5, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Tcell_SzaboModuleFeaturePlot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Tcell_SzaboModuleFeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("Szabo_Treg1","Szabo_CD4NV.CM1","Szabo_CD4.CD81",
                                                   "Szabo_IFNResponse1","Szabo_Proliferation1",
                                         "Szabo_CD8Cytotoxic1","Szabo_CD8Cytokine1")), pt.size = 0.1, y.max = 1)
ggsave("Tcell_Violin_Plot_res02_SzaboModules.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_SzaboModules.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("Cano_Gamez_Orange_Resting1","Cano_Gamez_Green_Resting1","Cano_Gamez_Orange1", "Cano_Gamez_Blue1","Cano_Gamez_Green1"), 
            pt.size = 0.5, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Tcell_CanoModuleFeaturePlot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Tcell_CanoModuleFeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("Cano_Gamez_Orange1","Cano_Gamez_Green1","Cano_Gamez_Green_High1","Cano_Gamez_Green_Low1")), pt.size = 0.1, y.max = 1.5)
ggsave("Tcell_Violin_Plot_res02_CanoModules.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_CanoModules.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Miao_CD4_naive <- list(gene.sets$Miao_CD4_naive)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_CD4_naive, name = "Miao_CD4_naive", seed.use = 1988)

Miao_CD8_naive <- list(gene.sets$Miao_CD8_naive)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_CD8_naive, name = "Miao_CD8_naive", seed.use = 1988)

Miao_Cytotoxic <- list(gene.sets$Miao_Cytotoxic)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Cytotoxic, name = "Miao_Cytotoxic", seed.use = 1988)

Miao_Exhausted <- list(gene.sets$Miao_Exhausted)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Exhausted, name = "Miao_Exhausted", seed.use = 1988)

Miao_Tr1 <- list(gene.sets$Miao_Tr1)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Tr1, name = "Miao_Tr1", seed.use = 1988)

Miao_nTreg <- list(gene.sets$Miao_nTreg)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_nTreg, name = "Miao_nTreg", seed.use = 1988)

Miao_iTreg <- list(gene.sets$Miao_iTreg)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_iTreg, name = "Miao_iTreg", seed.use = 1988)

Miao_Th1 <- list(gene.sets$Miao_Th1)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Th1, name = "Miao_Th1", seed.use = 1988)

Miao_Th2 <- list(gene.sets$Miao_Th2)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Th2, name = "Miao_Th2", seed.use = 1988)

Miao_Th17 <- list(gene.sets$Miao_Th17)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Th17, name = "Miao_Th17", seed.use = 1988)

Miao_Tfh <- list(gene.sets$Miao_Tfh)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Tfh, name = "Miao_Tfh", seed.use = 1988)

Miao_Central_memory <- list(gene.sets$Miao_Central_memory)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Central_memory, name = "Miao_Central_memory", seed.use = 1988)

Miao_Effector_memory <- list(gene.sets$Miao_Effector_memory)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Effector_memory, name = "Miao_Effector_memory", seed.use = 1988)

Miao_NKT <- list(gene.sets$Miao_NKT)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_NKT, name = "Miao_NKT", seed.use = 1988)

Miao_MAIT <- list(gene.sets$Miao_MAIT)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_MAIT, name = "Miao_MAIT", seed.use = 1988)

Miao_NK <- list(gene.sets$Miao_NK)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_NK, name = "Miao_NK", seed.use = 1988)

Miao_Gamma_delta <- list(gene.sets$Miao_Gamma_delta)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_Gamma_delta, name = "Miao_Gamma_delta", seed.use = 1988)

Miao_CD4_T <- list(gene.sets$Miao_CD4_T)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_CD4_T, name = "Miao_CD4_T", seed.use = 1988)

Miao_CD8_T <- list(gene.sets$Miao_CD8_T)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Miao_CD8_T, name = "Miao_CD8_T", seed.use = 1988)



head(seurat.all.data@meta.data)
unique(colnames(seurat.all.data@meta.data))

FeaturePlot(seurat.all.data, features = c("Miao_CD4_naive1","Miao_CD8_naive1","Miao_Cytotoxic1","Miao_Exhausted1",
                                          "Miao_Tr11","Miao_nTreg1","Miao_iTreg1","Miao_Th11","Miao_Th21","Miao_Th171"), pt.size = 0.5, 
            reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Tcell_MiaoModuleFeaturePlot1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Tcell_MiaoModuleFeaturePlot1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("Miao_Tfh1","Miao_Central_memory1","Miao_Effector_memory1",
                                          "Miao_NKT1","Miao_NK1","Miao_Gamma_delta1","Miao_CD4_T1","Miao_CD8_T1"), 
            pt.size = 0.5, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Tcell_MiaoModuleFeaturePlot2.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Tcell_MiaoModuleFeaturePlot2.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("Miao_CD4_naive1","Miao_CD8_naive1","Miao_Cytotoxic1","Miao_Exhausted1",
                                          "Miao_Tr11","Miao_nTreg1","Miao_iTreg1","Miao_Th11","Miao_Th21","Miao_Th171")), pt.size = 0.1, y.max = 1)
ggsave("Tcell_Violin_Plot_res02_MiaoModules1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_MiaoModules1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Miao_Tfh1","Miao_Central_memory1","Miao_Effector_memory1",
                                          "Miao_NKT1","Miao_NK1","Miao_Gamma_delta1","Miao_CD4_T1","Miao_CD8_T1")), pt.size = 0.1, y.max = 1)
ggsave("Tcell_Violin_Plot_res02_MiaoModules2.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_MiaoModules2.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

seurat.all.factors <- seurat.all.data@meta.data$ArrayName
plate.table = table(seurat.all.factors, as.vector(seurat.all.data@meta.data$RNA_snn_res.0.6))
head(plate.table)

write.table(plate.table,file="cellspercluster_res06.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp cellspercluster_res06.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_seurat_all_data_TcellFigure_correct.rds")
system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure_correct.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


metadata.forPCA <- seurat.all.data@meta.data
head(metadata.forPCA)

unique(metadata.forPCA$Replicate)

write.table(metadata.forPCA,file="modulescores_metadata.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp modulescores_metadata.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


unique(colnames(metadata.forPCA))
dim(metadata.forPCA)

metadata.forPCA.annotation <- metadata.forPCA[,1:54]

metadata.forPCA.annotation$percent.mt <- NULL
metadata.forPCA.annotation$seurat_clusters <- NULL
metadata.forPCA.annotation$pANN_0.25_0.005_1958 <- NULL
metadata.forPCA.annotation$DF.classifications_0.25_0.005_1958 <- NULL
metadata.forPCA.annotation$DF.classifications_0.25_0.005_1838 <- NULL
metadata.forPCA.annotation$original.clustering.1p2 <- NULL
metadata.forPCA.annotation$firstpassres2 <- NULL
metadata.forPCA.annotation$RNA_snn_res.2 <- NULL
metadata.forPCA.annotation$idents_res2 <- NULL
metadata.forPCA.annotation$RNA_snn_res.0.4 <- NULL
metadata.forPCA.annotation$RNA_snn_res.0.1 <- NULL

head(metadata.forPCA.annotation)


dim(metadata.forPCA)

module.scores <- metadata.forPCA[,25:55]

head(module.scores)
dim(module.scores)

row.annos <- cbind(metadata.forPCA.annotation$Hours,metadata.forPCA.annotation$Compartment)


row.annos <- as.data.frame(cbind(row.names(metadata.forPCA.annotation),row.annos))

row.names(row.annos) <- row.annos$V1
row.annos$V1 <- NULL
colnames(row.annos) <- c("Hours","Compartment")
head(row.annos)


row.annos <- cbind(row.annos,metadata.forPCA.annotation$RNA_snn_res.0.6)

#row.annos <- cbind(row.annos,metadata.forPCA.annotation$idents_res2_consensus)

row.annos <- cbind(row.annos, metadata.forPCA.annotation$Replicate)

colnames(row.annos) <- c("Hours","Compartment","Cluster", "Replicate")
head(row.annos)


unique(row.annos$Cluster) # 9 colors 
unique(row.annos$Hours) # 9 colors 
unique(row.annos$Compartment) # 9 colors 
unique(row.annos$Replicate) # 7 colors 

row.annos$Hours <- paste("Hours", row.annos$Hours, sep="_")

row.annos$Cluster <- paste("Cluster", row.annos$Cluster, sep="_")

row.annos$Replicate <- paste("Replicate", row.annos$Replicate, sep="_")

head(row.annos)

row.annos[order(row.annos$Compartment),]
row.annos[order(row.annos$Hours),]

dim(row.annos)

dim(module.scores)

unique(row.annos$Replicate)

cols.reps <- brewer.pal(7,"Paired")

cols.reps

colors_complex = list(Hours = c("Hours_050"="lightsteelblue1","Hours_066"="darkslategray1","Hours_073"="cadetblue2",
                                "Hours_089"="dodgerblue2","Hours_112"="blue1","Hours_137"="darkblue","Hours_500"="gold"),
                      Compartment = c("BLD"="brown","HEF"="darkgray"),
                      Cluster = c("Cluster_0" = "antiquewhite3","Cluster_1" = "indianred1","Cluster_2" = "hotpink2","Cluster_3" = "darkgoldenrod",
                    "Cluster_4" = "springgreen4","Cluster_5" = "darkolivegreen","Cluster_6" = "lightseagreen","Cluster_7" = "deepskyblue",
                    "Cluster_8" = "royalblue","Cluster_9" = "blue4","Cluster_10" = "purple","Cluster_11" = "darkorchid4"),
                     Replicate = c("Replicate_R1"= "orange", "Replicate_R2" = "orange", "Replicate_00" = "gold", "Replicate_01" = "gold",
                                   "Replicate_02" = "gold", "Replicate_03" = "gold", "Replicate_PT" = "orange"))


top_ha = HeatmapAnnotation(df = row.annos, col = colors_complex)


top_ha

Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE,
        cluster_columns = FALSE, show_row_names = TRUE)

Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,
        cluster_columns = FALSE, show_row_names = TRUE)

HA <- Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,
        cluster_columns = FALSE, show_row_names = TRUE) 
pdf("Tcell_ModuleHeatmap_stimannos.pdf")
print(HA)
dev.off()
system(paste0("gsutil cp Tcell_ModuleHeatmap_stimannos.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


HA_nameannos <- Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,
        cluster_columns = FALSE, show_row_names = TRUE)
pdf("Tcell_ModuleHeatmap_nameannos.pdf")
print(HA_nameannos)
dev.off()
system(paste0("gsutil cp Tcell_ModuleHeatmap_nameannos.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_seurat_all_data_TcellFigure.rds")
system(paste0("gsutil cp seuratobj_seurat_all_data_TcellFigure.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


all.pca <- prcomp(module.scores)


ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, var.axes = FALSE,varname.adjust = 2) + scale_color_manual(values=colors.cluster.new)
ggsave("TcellBiplot_nolabels.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp TcellBiplot_nolabels.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2) + scale_color_manual(values=colors.cluster.new)
ggsave("TcellBiplot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp TcellBiplot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2,var.axes = FALSE, choices = 2:3) + scale_color_manual(values=colors.cluster.new)
ggsave("TcellBiplot_PCs2and3.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp TcellBiplot_PCs2and3.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2,var.axes = FALSE, choices = 3:4) + scale_color_manual(values=colors.cluster.new)
ggsave("TcellBiplot_PCs3and4.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp TcellBiplot_PCs3and4.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


var.loadings <- all.pca$rotation


write.table(var.loadings,file="variable_loadings.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp variable_loadings.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


module.loading.PC1 <- as.data.frame(all.pca$rotation[,1])
colnames(module.loading.PC1) <- "PC1"
module.loading.PC1$modules <- row.names(module.loading.PC1)
module.loading.PC1.ordered <- arrange(module.loading.PC1,module.loading.PC1$PC1)
head(module.loading.PC1.ordered)
tail(module.loading.PC1.ordered)

module.loading.PC2 <- as.data.frame(all.pca$rotation[,2])
colnames(module.loading.PC2) <- "PC2"
module.loading.PC2$modules <- row.names(module.loading.PC2)
module.loading.PC2.ordered <- arrange(module.loading.PC2,module.loading.PC2$PC2)
head(module.loading.PC2.ordered)
tail(module.loading.PC2.ordered)

FeaturePlot(seurat.all.data, features = c("Szabo_CD4NV.CM1"), pt.size = 0.5, 
            reduction = "tsne", cols = c("white","black"))

FeaturePlot(seurat.all.data, features = c("Cano_Gamez_Orange_Resting1"), pt.size = 0.5, 
            reduction = "tsne", cols = c("blue", "yellow"))

# make feature plots 

FeaturePlot(seurat.all.data, features = c("Cano_Gamez_Orange_Resting1"), pt.size = 1, 
            reduction = "tsne", cols = c("yellow","blue"))
ggsave("ModuleFeaturePlot_Cano_Gamez_Orange_Resting1.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp ModuleFeaturePlot_Cano_Gamez_Orange_Resting1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("Cano_Gamez_Green_Resting1"), pt.size = 1, 
            reduction = "tsne", cols = c("yellow","blue"))
ggsave("ModuleFeaturePlot_Cano_Gamez_Green_Resting1.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp ModuleFeaturePlot_Cano_Gamez_Green_Resting1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("Szabo_CD4NV.CM1"), pt.size = 1, 
            reduction = "tsne", cols = c("yellow","blue"))
ggsave("ModuleFeaturePlot_Szabo_CD4NVCM1.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp ModuleFeaturePlot_Szabo_CD4NVCM1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


head(var.loadings)

## make violin plots of modules: Cano_Gamez_Green1, Miao_Cytotoxic1, Miao_iTreg1, Miao_CD4_naive1 (PC1), Cano_Gamez_Orange1, Miao_NK1 (PC2)

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("Cano_Gamez_Green_Resting1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Cano_Gamez_Green_Resting1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Cano_Gamez_Green_Resting1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Miao_Tr11")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Miao_Tr11.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Miao_Tr11.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Miao_CD4_naive1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Miao_CD4_naive1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Miao_CD4_naive1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Cano_Gamez_Orange_Resting1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Cano_Gamez_Orange_Resting1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Cano_Gamez_Orange_Resting1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Cano_Gamez_Green1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Cano_Gamez_Green1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Cano_Gamez_Green1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("Szabo_CD4NV.CM1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_Szabo_CD4NVCM1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_Szabo_CD4NVCM1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



pcs.matrix = all.pca$x[,c(1:2)]


plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3)


rownames(as.data.frame(all.pca$rotation[,1]))

pdf("AllModule_PCA.pdf")
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3)
dev.off()
system(paste0("gsutil cp AllModule_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


var.loadings[2,][1]

plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1

#Cano_Gamez_Green_Resting1

pdf("Cano_Gamez_Green_Resting1_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
dev.off()
system(paste0("gsutil cp Cano_Gamez_Green_Resting1_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


var.loadings[7,][1]

#Szabo_CD4NV.CM1

pdf("Szabo_CD4NVCM1_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Szabo_CD4NVCM1_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


var.loadings[1,][2]

#Cano_Gamez_Orange_Resting1

pdf("Cano_Gamez_Orange_Resting1_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
dev.off()
system(paste0("gsutil cp Cano_Gamez_Orange_Resting1_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


dim(var.loadings)

pdf("AllmoduleArrows_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.3, col = "white")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[3,][1],var.loadings[3,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[4,][1],var.loadings[4,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[5,][1],var.loadings[5,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[6,][1],var.loadings[6,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
arrows(0,0,var.loadings[8,][1],var.loadings[8,][2],col = "grey", code = 0) 
arrows(0,0,var.loadings[9,][1],var.loadings[9,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[12,][1],var.loadings[12,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[13,][1],var.loadings[13,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[14,][1],var.loadings[14,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[15,][1],var.loadings[15,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[16,][1],var.loadings[16,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[17,][1],var.loadings[17,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[18,][1],var.loadings[18,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[19,][1],var.loadings[19,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[20,][1],var.loadings[20,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[21,][1],var.loadings[21,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[22,][1],var.loadings[22,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[23,][1],var.loadings[23,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[24,][1],var.loadings[24,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[25,][1],var.loadings[25,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[26,][1],var.loadings[26,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[27,][1],var.loadings[27,][2],col = "grey", code = 0) 
arrows(0,0,var.loadings[28,][1],var.loadings[28,][2],col = "grey", code = 0) 
arrows(0,0,var.loadings[29,][1],var.loadings[26,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[30,][1],var.loadings[27,][2],col = "grey", code = 0) 
arrows(0,0,var.loadings[31,][1],var.loadings[28,][2],col = "grey", code = 0) 

dev.off()
system(paste0("gsutil cp AllmoduleArrows_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)



head(row.annos)
unique(row.annos$Hours)

pdf("Patient.pdf", useDingbats = F)
cells.use.replicate = rownames(row.annos)[which(row.annos$Replicate=="Replicate_R1" | row.annos$Replicate=="Replicate_R2" | row.annos$Replicate=="Replicate_PT")]
plot(pcs.matrix[cells.use.replicate,1], pcs.matrix[cells.use.replicate,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="Patient", col="orange")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Patient.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)




pdf("Controls.pdf", useDingbats = F)
cells.use.replicate = rownames(row.annos)[which(row.annos$Replicate=="Replicate_00" | row.annos$Replicate=="Replicate_01" | row.annos$Replicate=="Replicate_02" | row.annos$Replicate=="Replicate_03")]
plot(pcs.matrix[cells.use.replicate,1], pcs.matrix[cells.use.replicate,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="Controls", col="gold")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Controls.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)




pdf("Cluster0_PCA.pdf", useDingbats = F)
cells.use.c0 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_0")]
plot(pcs.matrix[cells.use.c0,1], pcs.matrix[cells.use.c0,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 0", col="antiquewhite3")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster0_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster1_PCA.pdf", useDingbats = F)
cells.use.c1 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_1")]
plot(pcs.matrix[cells.use.c1,1], pcs.matrix[cells.use.c1,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 1", col="indianred1")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster1_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster2_PCA.pdf", useDingbats = F)
cells.use.c2 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_2")]
plot(pcs.matrix[cells.use.c2,1], pcs.matrix[cells.use.c2,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 2", col="hotpink2")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster2_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster3_PCA.pdf", useDingbats = F)
cells.use.c3 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_3")]
plot(pcs.matrix[cells.use.c3,1], pcs.matrix[cells.use.c3,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 3", col="darkgoldenrod")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster3_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster4_PCA.pdf", useDingbats = F)
cells.use.c4 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_4")]
plot(pcs.matrix[cells.use.c4,1], pcs.matrix[cells.use.c4,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 4", col="springgreen4")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster4_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster5_PCA.pdf", useDingbats = F)
cells.use.c5 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_5")]
plot(pcs.matrix[cells.use.c5,1], pcs.matrix[cells.use.c5,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 5", col="darkolivegreen")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster5_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster6_PCA.pdf", useDingbats = F)
cells.use.c6 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_6")]
plot(pcs.matrix[cells.use.c6,1], pcs.matrix[cells.use.c6,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 6", col="lightseagreen")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster6_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster7_PCA.pdf", useDingbats = F)
cells.use.c7 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_7")]
plot(pcs.matrix[cells.use.c7,1], pcs.matrix[cells.use.c7,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 7", col="deepskyblue")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster7_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster8_PCA.pdf", useDingbats = F)
cells.use.c8 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_8")]
plot(pcs.matrix[cells.use.c8,1], pcs.matrix[cells.use.c8,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 8", col="royalblue")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster8_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster9_PCA.pdf", useDingbats = F)
cells.use.c9 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_9")]
plot(pcs.matrix[cells.use.c9,1], pcs.matrix[cells.use.c9,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 9", col="blue4")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster9_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster10_PCA.pdf", useDingbats = F)
cells.use.c10 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_10")]
plot(pcs.matrix[cells.use.c10,1], pcs.matrix[cells.use.c10,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 10", col="purple")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster10_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


pdf("Cluster11_PCA.pdf", useDingbats = F)
cells.use.c11 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_11")]
plot(pcs.matrix[cells.use.c11,1], pcs.matrix[cells.use.c11,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="cluster 11", col="darkorchid4")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp Cluster11_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# By time

pdf("hematoma_50hrs_PCA.pdf", useDingbats = F)
cells.use.50hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_050")]
plot(pcs.matrix[cells.use.50hrs.hematoma,1], pcs.matrix[cells.use.50hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="50 hours, hematoma", col="lightsteelblue1")
dev.off()
system(paste0("gsutil cp hematoma_50hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("hematoma_66hrs_PCA.pdf", useDingbats = F)
cells.use.66hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_066")]
plot(pcs.matrix[cells.use.66hrs.hematoma,1], pcs.matrix[cells.use.66hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="66 hours, hematoma", col="darkslategray1")
dev.off()
system(paste0("gsutil cp hematoma_66hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("hematoma_73hrs_PCA.pdf", useDingbats = F)
cells.use.73hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_073")]
plot(pcs.matrix[cells.use.73hrs.hematoma,1], pcs.matrix[cells.use.73hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="73 hours, hematoma", col="cadetblue2")
dev.off()
system(paste0("gsutil cp hematoma_73hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("hematoma_89hrs_PCA.pdf", useDingbats = F)
cells.use.89hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_089")]
plot(pcs.matrix[cells.use.89hrs.hematoma,1], pcs.matrix[cells.use.89hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="89 hours, hematoma", col="dodgerblue2")
dev.off()
system(paste0("gsutil cp hematoma_89hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("hematoma_112hrs_PCA.pdf", useDingbats = F)
cells.use.112hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_112")]
plot(pcs.matrix[cells.use.112hrs.hematoma,1], pcs.matrix[cells.use.112hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="112 hours, hematoma", col="blue1")
dev.off()
system(paste0("gsutil cp hematoma_112hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("hematoma_137hrs_PCA.pdf", useDingbats = F)
cells.use.137hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_137")]
plot(pcs.matrix[cells.use.137hrs.hematoma,1], pcs.matrix[cells.use.137hrs.hematoma,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="137 hours, hematoma", col="darkblue")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp hematoma_137hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_50hrs_PCA.pdf", useDingbats = F)
cells.use.50hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_050")]
plot(pcs.matrix[cells.use.50hrs.blood,1], pcs.matrix[cells.use.50hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="50 hours, blood", col="lightsteelblue1")
dev.off()
system(paste0("gsutil cp blood_50hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_66hrs_PCA.pdf", useDingbats = F)
cells.use.66hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_066")]
plot(pcs.matrix[cells.use.66hrs.blood,1], pcs.matrix[cells.use.66hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="66 hours, blood", col="darkslategray1")
dev.off()
system(paste0("gsutil cp blood_66hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_73hrs_PCA.pdf", useDingbats = F)
cells.use.73hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_073")]
plot(pcs.matrix[cells.use.73hrs.blood,1], pcs.matrix[cells.use.73hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="73 hours, blood", col="cadetblue2")
dev.off()
system(paste0("gsutil cp blood_73hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_89hrs_PCA.pdf", useDingbats = F)
cells.use.89hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_089")]
plot(pcs.matrix[cells.use.89hrs.blood,1], pcs.matrix[cells.use.89hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="89 hours, blood", col="dodgerblue2")
dev.off()
system(paste0("gsutil cp blood_89hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_112hrs_PCA.pdf", useDingbats = F)
cells.use.112hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_112")]
plot(pcs.matrix[cells.use.112hrs.blood,1], pcs.matrix[cells.use.112hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="112 hours, blood", col="blue1")
dev.off()
system(paste0("gsutil cp blood_112hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_137hrs_PCA.pdf", useDingbats = F)
cells.use.137hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_137")]
plot(pcs.matrix[cells.use.137hrs.blood,1], pcs.matrix[cells.use.137hrs.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="137 hours, blood", col="darkblue")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp blood_137hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)

pdf("blood_ControlsandPT_PCA.pdf", useDingbats = F)
cells.use.control.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_500")]
plot(pcs.matrix[cells.use.control.blood,1], pcs.matrix[cells.use.control.blood,2], xlim=c(-1,2.1), ylim=c(-2,1.5), pch=19, cex=.7, main="controls+PT, blood", col="gold")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2]) #Cano_Gamez_Orange_Resting1
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2]) #Cano_Gamez_Green_Resting1
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2]) #Szabo_CD4NVCM1
dev.off()
system(paste0("gsutil cp blood_ControlsandPT_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# NOTE: Need to clean up effectorness calculations (exported and merged outside Terra for now)

#library(monocle)
library(monocle3)


#Tcells.monocle <- importCDS(seurat.all.data) ## function still doesn't work apparently, code below accomplishes the same thing 


data <- as(as.matrix(seurat.all.data@assays$RNA@data), 'sparseMatrix')


pd <- new('AnnotatedDataFrame', data = seurat.all.data@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)


monocle_cds <- new_cell_data_set(data,
                         cell_metadata = seurat.all.data@meta.data,
                         gene_metadata = fData)

monocle_cds <- preprocess_cds(monocle_cds, num_dim = 100)

monocle_cds <- reduce_dimension(monocle_cds, reduction_method = 'UMAP')


monocle_cds <- cluster_cells(monocle_cds, reduction_method = 'UMAP')


monocle_cds <- learn_graph(monocle_cds)


#monocle_cds <- order_cells(monocle_cds)


plot_cells(monocle_cds, color_cells_by = "partition")
ggsave("Monocle3_partitionplot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Monocle3_partitionplot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


plot_cells(monocle_cds, color_cells_by = "RNA_snn_res.0.2", label_leaves=TRUE,
           label_branch_points=TRUE)
ggsave("Monocle3_NTclusterID.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Monocle3_NTclusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


head(colData(monocle_cds))

# Note: this function doesn't work 
get_earliest_principal_node <- function(cds, time_bin="RNA_snn_res.0.2"){
  cell_ids <- which(colData(cds)[, "6"] == time_bin)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

clusters(monocle_cds)

# run pseudotime over all cells
root_group = colnames(monocle_cds)[clusters(monocle_cds) == 10]


monocle_cds <- order_cells(monocle_cds, root_cells=root_group)

# Root 9 looks like it starts opposite of effector scores, near clustr 3 in seurat (use 9 for partition 1)
# Root 7 looks like it has lower pseudotime on NK cluster at top, need reversed for parition 2 --> 10! 
plot_cells(monocle_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
#ggsave("Monocle3_pseudotime_root9_final.pdf", useDingbats = F, height = 12, width = 16, units = "in")
#system(paste0("gsutil cp Monocle3_pseudotime_root9_final.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#extract partition 1
partion2 <- as.data.frame(pseudotime(monocle_cds))

tail(partion2)

#ICH.EffCalc.Part1 <- (pseudotime(monocle_cds) - min(pseudotime(monocle_cds)))/(max(pseudotime(monocle_cds))-min(pseudotime(monocle_cds)))

#ICH.Eff <- (pData(monocle_cds)$Pseudotime - min(pData(monocle_cds)$Pseudotime))/(max(pData(monocle_cds)$pseudotime)-min(pData(monocle_cds)$pseudotime))


tail(ICH.EffCalc.Part1)

colnames(partion2) <- c("Pseudotime_partition2")


write.table(partion2,file="ICH_eff_partition2.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp ICH_eff_partition2.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)




tail(paritiion1)

partition2 <- pseudotime(monocle_cds)


tail(partition2)

ICH.EffCalc.Part2 <- as.data.frame((pseudotime(monocle_cds) - min(pseudotime(monocle_cds)))/(max(pseudotime(monocle_cds))-min(pseudotime(monocle_cds))))


dim(ICH.EffCalc.Part2)

colnames(ICH.EffCalc.Part2) <- c("Pseudotime_partition2")
tail(ICH.EffCalc.Part2)

write.table(ICH.EffCalc.Part2,file="ICH_eff_partition2.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp ICH_eff_partition2.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# note: before calculating effectorness, need to remove the cells that are not included?
#ICH.Eff <- (pData(monocle_cds)$Pseudotime - min(pData(monocle_cds)$Pseudotime))/(max(pData(monocle_cds)$pseudotime)-min(pData(monocle_cds)$pseudotime))



pData(monocle_cds)

head(ICH.Eff)

help(export_cds)

help(order_cells)

effectorness <- read.csv(pipe('gsutil cat gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/merged_psuedotime.csv'))


head(effectorness)

rownames(effectorness) <- effectorness$X
effectorness$X <- NULL
effectorness$Order <- NULL

seurat.all.data <- AddMetaData(seurat.all.data, effectorness)

head(seurat.all.data@meta.data)

FeaturePlot(seurat.all.data, features = c("effectorness_merged"), 
            pt.size = 1.4, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Effectorness.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Effectorness.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c(features = c("effectorness_merged")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res02_Effectorness.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_Effectorness.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("pseudotime_merged"), 
            pt.size = 0.5, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_pseudotime.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_pseudotime.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c(features = c("pseudotime_merged")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res02_Pseudotime.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res02_Pseudotime.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


colnames(seurat.all.data@meta.data)
unique(seurat.all.data@meta.data$idents_res2_consensus)

Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.nonCD8s <- subset(seurat.all.data, idents = c("activated T cells",
"proliferating T cells",
"heat shock response T cells",
"T cells 4",
"T cells 2",
"ribosomal-enriched T cells",
"XIST+ T cells",
"T cells 1",
"T cells 5",
"T cells 6",
"T cells 7"))


unique(suerat.all.data.nonCD8s@meta.data$Compartment)

Idents(object=suerat.all.data.nonCD8s) <- "Hours"
suerat.all.data.nonCD8s.050 <- subset(suerat.all.data.nonCD8s, idents = c("050"))
suerat.all.data.nonCD8s.066 <- subset(suerat.all.data.nonCD8s, idents = c("066"))
suerat.all.data.nonCD8s.073 <- subset(suerat.all.data.nonCD8s, idents = c("073"))
suerat.all.data.nonCD8s.089 <- subset(suerat.all.data.nonCD8s, idents = c("089"))
suerat.all.data.nonCD8s.112 <- subset(suerat.all.data.nonCD8s, idents = c("112"))
suerat.all.data.nonCD8s.137 <- subset(suerat.all.data.nonCD8s, idents = c("137"))

Idents(object=suerat.all.data.nonCD8s.050) <- "Compartment"
suerat.all.data.nonCD8s.050.BLD <- subset(suerat.all.data.nonCD8s.050, idents = c("BLD"))
suerat.all.data.nonCD8s.050.HEF <- subset(suerat.all.data.nonCD8s.050, idents = c("HEF"))

Idents(object=suerat.all.data.nonCD8s.066) <- "Compartment"
suerat.all.data.nonCD8s.066.BLD <- subset(suerat.all.data.nonCD8s.066, idents = c("BLD"))
suerat.all.data.nonCD8s.066.HEF <- subset(suerat.all.data.nonCD8s.066, idents = c("HEF"))

Idents(object=suerat.all.data.nonCD8s.073) <- "Compartment"
suerat.all.data.nonCD8s.073.BLD <- subset(suerat.all.data.nonCD8s.073, idents = c("BLD"))
suerat.all.data.nonCD8s.073.HEF <- subset(suerat.all.data.nonCD8s.073, idents = c("HEF"))

Idents(object=suerat.all.data.nonCD8s.089) <- "Compartment"
suerat.all.data.nonCD8s.089.BLD <- subset(suerat.all.data.nonCD8s.089, idents = c("BLD"))
suerat.all.data.nonCD8s.089.HEF <- subset(suerat.all.data.nonCD8s.089, idents = c("HEF"))

Idents(object=suerat.all.data.nonCD8s.112) <- "Compartment"
suerat.all.data.nonCD8s.112.BLD <- subset(suerat.all.data.nonCD8s.112, idents = c("BLD"))
suerat.all.data.nonCD8s.112.HEF <- subset(suerat.all.data.nonCD8s.112, idents = c("HEF"))

Idents(object=suerat.all.data.nonCD8s.137) <- "Compartment"
suerat.all.data.nonCD8s.137.BLD <- subset(suerat.all.data.nonCD8s.137, idents = c("BLD"))
suerat.all.data.nonCD8s.137.HEF <- subset(suerat.all.data.nonCD8s.137, idents = c("HEF"))

# 50 
count.data.050.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.050.BLD, slot = "counts"))
count.data.050.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.050.HEF, slot = "counts"))

count.data.050.sum <- cbind(rowSums(count.data.050.BLD),rowSums(count.data.050.HEF))

write.table(count.data.050.sum,file="Pseudobulk_050hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_050hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


# 66
count.data.066.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.066.BLD, slot = "counts"))
count.data.066.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.066.HEF, slot = "counts"))

count.data.066.sum <- cbind(rowSums(count.data.066.BLD),rowSums(count.data.066.HEF))

write.table(count.data.066.sum,file="Pseudobulk_066hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_066hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#
count.data.073.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.073.BLD, slot = "counts"))
count.data.073.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.073.HEF, slot = "counts"))

count.data.073.sum <- cbind(rowSums(count.data.073.BLD),rowSums(count.data.073.HEF))

write.table(count.data.073.sum,file="Pseudobulk_073hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_073hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#
count.data.089.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.089.BLD, slot = "counts"))
count.data.089.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.089.HEF, slot = "counts"))

count.data.089.sum <- cbind(rowSums(count.data.089.BLD),rowSums(count.data.089.HEF))

write.table(count.data.089.sum,file="Pseudobulk_089hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_089hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#
count.data.112.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.112.BLD, slot = "counts"))
count.data.112.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.112.HEF, slot = "counts"))

count.data.112.sum <- cbind(rowSums(count.data.112.BLD),rowSums(count.data.112.HEF))

write.table(count.data.112.sum,file="Pseudobulk_112hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_112hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


#
count.data.137.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.137.BLD, slot = "counts"))
count.data.137.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.nonCD8s.137.HEF, slot = "counts"))

count.data.137.sum <- cbind(rowSums(count.data.137.BLD),rowSums(count.data.137.HEF))

write.table(count.data.137.sum,file="Pseudobulk_137hrs_Tcells.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_137hrs_Tcells.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


colnames(seurat.all.data@meta.data)

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster0 <- subset(seurat.all.data, idents = c("0"))

suerat.all.data.subset.cluster0

unique(suerat.all.data.subset.cluster0@meta.data$RNA_snn_res.0.6)

length(suerat.all.data.subset.cluster0@meta.data$Compartment == "BLD")

length(suerat.all.data.subset.cluster0@meta.data$Compartment == "HEF")

seurat.factors.cluster1 <- suerat.all.data.subset.cluster0@meta.data$RNA_snn_res.0.6
plate.table.cluster1 = table(seurat.factors.cluster1, as.vector(suerat.all.data.subset.cluster0@meta.data$Compartment))
head(plate.table.cluster1)

Idents(object=suerat.all.data.subset.cluster0) <- "Compartment"

seurat.all.data.markers.cluster0 <- FindAllMarkers(suerat.all.data.subset.cluster0, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 230, random.seed = 1988)

write.table(seurat.all.data.markers.cluster0,file="seurat.all.data.markers.cluster0.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster0.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


head(seurat.all.data.markers.cluster0)

top10.cluster0 <- seurat.all.data.markers.cluster0 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster0,features=top10.cluster0$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster0.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster0.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster2 <- subset(seurat.all.data, idents = c("2"))


seurat.factors.cluster2 <- suerat.all.data.subset.cluster2@meta.data$RNA_snn_res.0.6
plate.table.cluster2 = table(seurat.factors.cluster2, as.vector(suerat.all.data.subset.cluster2@meta.data$Compartment))
head(plate.table.cluster2)

Idents(object=suerat.all.data.subset.cluster2) <- "Compartment"

seurat.all.data.markers.cluster2 <- FindAllMarkers(suerat.all.data.subset.cluster2, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 412, random.seed = 1988)

write.table(seurat.all.data.markers.cluster2,file="seurat.all.data.markers.cluster2.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster2.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster2 <- seurat.all.data.markers.cluster2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster2,features=top10.cluster2$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster2.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster2.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster3 <- subset(seurat.all.data, idents = c("3"))


seurat.factors.cluster3 <- suerat.all.data.subset.cluster3@meta.data$RNA_snn_res.0.6
plate.table.cluster3 = table(seurat.factors.cluster3, as.vector(suerat.all.data.subset.cluster3@meta.data$Compartment))
head(plate.table.cluster3)

Idents(object=suerat.all.data.subset.cluster3) <- "Compartment"

seurat.all.data.markers.cluster3 <- FindAllMarkers(suerat.all.data.subset.cluster3, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 385, random.seed = 1988)

write.table(seurat.all.data.markers.cluster3,file="seurat.all.data.markers.cluster3.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster3.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster3 <- seurat.all.data.markers.cluster3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster3,features=top10.cluster3$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster3.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster3.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster4 <- subset(seurat.all.data, idents = c("4"))


seurat.factors.cluster4 <- suerat.all.data.subset.cluster4@meta.data$RNA_snn_res.0.6
plate.table.cluster4 = table(seurat.factors.cluster4, as.vector(suerat.all.data.subset.cluster4@meta.data$Compartment))
head(plate.table.cluster4)


Idents(object=suerat.all.data.subset.cluster4) <- "Compartment"

seurat.all.data.markers.cluster4 <- FindAllMarkers(suerat.all.data.subset.cluster4, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 231, random.seed = 1988)

write.table(seurat.all.data.markers.cluster4,file="seurat.all.data.markers.cluster4.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster4.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster4 <- seurat.all.data.markers.cluster4 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster4,features=top10.cluster4$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster4.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster4.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster7 <- subset(seurat.all.data, idents = c("7"))


seurat.factors.cluster7 <- suerat.all.data.subset.cluster7@meta.data$RNA_snn_res.0.6
plate.table.cluster7 = table(seurat.factors.cluster7, as.vector(suerat.all.data.subset.cluster7@meta.data$Compartment))
plate.table.cluster7

Idents(object=suerat.all.data.subset.cluster7) <- "Compartment"

seurat.all.data.markers.cluster7 <- FindAllMarkers(suerat.all.data.subset.cluster7, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 511, random.seed = 1988)

write.table(seurat.all.data.markers.cluster7,file="seurat.all.data.markers.cluster7.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster7.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster7 <- seurat.all.data.markers.cluster7 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster7,features=top10.cluster7$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster7.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster7.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster9 <- subset(seurat.all.data, idents = c("9"))


seurat.factors.cluster9 <- suerat.all.data.subset.cluster9@meta.data$RNA_snn_res.0.6
plate.table.cluster9 = table(seurat.factors.cluster9, as.vector(suerat.all.data.subset.cluster9@meta.data$Compartment))
plate.table.cluster9

Idents(object=suerat.all.data.subset.cluster9) <- "Compartment"

seurat.all.data.markers.cluster9 <- FindAllMarkers(suerat.all.data.subset.cluster9, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 211, random.seed = 1988)

write.table(seurat.all.data.markers.cluster9,file="seurat.all.data.markers.cluster9.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster9.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster9 <- seurat.all.data.markers.cluster9 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster9,features=top10.cluster9$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster9.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster9.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster10 <- subset(seurat.all.data, idents = c("10"))


seurat.factors.cluster10 <- suerat.all.data.subset.cluster10@meta.data$RNA_snn_res.0.6
plate.table.cluster10 = table(seurat.factors.cluster10, as.vector(suerat.all.data.subset.cluster10@meta.data$Compartment))
plate.table.cluster10


Idents(object=suerat.all.data.subset.cluster10) <- "Compartment"

seurat.all.data.markers.cluster10 <- FindAllMarkers(suerat.all.data.subset.cluster10, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 29, random.seed = 1988)

write.table(seurat.all.data.markers.cluster10,file="seurat.all.data.markers.cluster10.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster10.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)


top10.cluster10 <- seurat.all.data.markers.cluster10 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster10,features=top10.cluster10$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster10.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster10.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200403_TcellAnalysis/"), intern=TRUE)















