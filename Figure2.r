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

library("gplots")

project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
project
workspace
bucket

list.files(path = ".")

# use below path for first run through analysis, the second for subsequent ones
# RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/*.rds"

RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/*.rds"
system(sprintf("gsutil ls %s", RDS_Files), intern=T)

system("mkdir $(pwd)/RDS_Files")
system(sprintf("gsutil -m cp %s $(pwd)/RDS_Files/", RDS_Files))

list.files(path = "./RDS_Files")

#seurat.all.data <- readRDS("./RDS_Files/seuratobj_alldata_filtered_myeloid.rds") #import first
seurat.all.data <- readRDS("./RDS_Files/seuratobj_alldata_filtered_monos.rds") #import after removing grans


seurat.all.data

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
system(paste0("gsutil cp Variable_Genes_Merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# scale and PCA, can regress stuff out here too (do not regress for now)
seurat.all.data <- ScaleData(seurat.all.data, features = rownames(seurat.all.data))
seurat.all.data <- RunPCA(seurat.all.data, features = VariableFeatures(object = seurat.all.data))

ElbowPlot(seurat.all.data)
ggsave("Elbow_merge.pdf", useDingbats = F)
system(paste0("gsutil cp Elbow_merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


DimHeatmap(seurat.all.data, dims = 1:20, cells = 500, balanced = TRUE)
ggsave("PC_heatmaps.pdf", useDingbats = F, height = 10, width = 10, units = "in")
system(paste0("gsutil cp PC_heatmaps.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


head(seurat.all.data@meta.data)

Idents(object=seurat.all.data) <- "ArrayName"

DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_array.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_array.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "Platform"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_platform.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_platform.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# PCA plots colored by compartment
Idents(object=seurat.all.data) <- "Compartment"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# PCA plots colored by SNN2 cluster 
Idents(object=seurat.all.data) <- "idents_res2_consensus"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_consensus.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_consensus.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_consensus.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_consensus.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# PCA plots colored by time
Idents(object=seurat.all.data) <- "Hours"
DimPlot(seurat.all.data, reduction = "pca", dims = c(1,2))
ggsave("PC1_2_hours.pdf", useDingbats = F)
system(paste0("gsutil cp PC1_2_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "pca", dims = c(3,4))
ggsave("PC3_4_hours.pdf", useDingbats = F)
system(paste0("gsutil cp PC3_4_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# remove granulocytes, re-import and remake figures above
#Idents(object=seurat.all.data) <- "idents_res2_consensus"
#suerat.all.data.subset <- subset(seurat.all.data, idents = c("monocyte/macrophage 1",
#"CD16 non-classical monocytes",
#"CD1C+ DCs",
#"monocyte/macrophage 6",
#"monocyte/macrophage 2",
#"macrophage 4",
#"macrophage 1",
#"monocyte/macrophage 3",
#"monocyte/macrophage 4",
#"monocyte-RBCs",
#"monocyte/macrophage 5",
#"macrophage 3",
#"macrophage 2"))
#saveRDS(suerat.all.data.subset,file="seuratobj_alldata_filtered_monos.rds")
#system(paste0("gsutil cp seuratobj_alldata_filtered_monos.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



## Jackstraw sig PCs (takes a really long time to run - all 20 PCs are sig!)
#seurat.all.data <- JackStraw(seurat.all.data, num.replicate = 100)
#seurat.all.data <- ScoreJackStraw(seurat.all.data, dims = 1:20)
#JackStrawPlot(seurat.all.data, dims = 1:20)
#ggsave("JackstrawPlot.pdf", useDingbats = F)
#system(paste0("gsutil cp JackstrawPlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


seurat.all.data <- RunTSNE(seurat.all.data, reduction = "pca", dims = c(1:20))
seurat.all.data <- FindNeighbors(seurat.all.data, reduction = "pca", dims = c(1:20))


# QC feature plots
FeaturePlot(seurat.all.data, features = c("percent.mt"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_percentmito.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_percentmito.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# colored by other features
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform)
ggsave("UMAP_platform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform")
ggsave("UMAP_platform_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment)
ggsave("UMAP_compartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment")
ggsave("UMAP_compartment_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8)
ggsave("UMAP_arrayname.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours)
ggsave("UMAP_hours.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours")
ggsave("UMAP_hours_splitbyhrs.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate")
ggsave("UMAP_hours_splitbyreplicate.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform")
ggsave("UMAP_hours_splitbyplatform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment")
ggsave("UMAP_hours_splitbycompartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# cell ID
DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8)
ggsave("UMAP_cellID.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Hours")
ggsave("UMAP_cellID_splitbyhrs.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyhrs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Replicate")
ggsave("UMAP_cellID_splitbyreplicate.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Platform")
ggsave("UMAP_cellID_splitbyplatform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Compartment")
ggsave("UMAP_cellID_splitbycompartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbycompartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



# colored by other features (no legend)
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform) + NoLegend()
ggsave("UMAP_platform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform") + NoLegend()
ggsave("UMAP_platform_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment) + NoLegend()
ggsave("UMAP_compartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment") + NoLegend()
ggsave("UMAP_compartment_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8) + NoLegend()
ggsave("UMAP_arrayname_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours) + NoLegend()
ggsave("UMAP_hours_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours") + NoLegend()
ggsave("UMAP_hours_splitbyhrs_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate") + NoLegend()
ggsave("UMAP_hours_splitbyreplicate_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform") + NoLegend()
ggsave("UMAP_hours_splitbyplatform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment") + NoLegend()
ggsave("UMAP_hours_splitbycompartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# cell ID
DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8)+ NoLegend()
ggsave("UMAP_cellID_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Hours")+ NoLegend()
ggsave("UMAP_cellID_splitbyhrs_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyhrs_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Replicate")+ NoLegend()
ggsave("UMAP_cellID_splitbyreplicate_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Platform")+ NoLegend()
ggsave("UMAP_cellID_splitbyplatform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "idents_res2_consensus", pt.size = 0.8, split.by = "Compartment")+ NoLegend()
ggsave("UMAP_cellID_splitbycompartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_cellID_splitbycompartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



# color tSNE with 

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
system(paste0("gsutil cp tsne_time_collapse_PThighlight.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.array.pt.ctl) + NoLegend()
ggsave("tsne_time_collapse_PThighlight_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_time_collapse_PThighlight_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)




# cluster IDs (0.4, 0.3, 0.5, 0.1, 0.2)
seurat.all.data <- FindClusters(seurat.all.data, resolution = 0.6)


# all together
DimPlot(seurat.all.data, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colors.cluster.new)  
ggsave("tsne_clusterID_res06.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = FALSE, pt.size = 0.5, cols = colors.cluster.new) + NoLegend()
ggsave("tsne_clusterID__res06_no_legend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID__res06_no_legend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#split by compartment
DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment", cols = colors.cluster.new)  
ggsave("tsne_clusterID_res06_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment", cols = colors.cluster.new) + NoLegend()
ggsave("tsne_clusterID_res06_no_legend_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_res06_no_legend_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



head(seurat.all.data@meta.data)

## Try 0.6
Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
seurat.all.data.markers <- FindAllMarkers(seurat.all.data, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5)

write.table(seurat.all.data.markers,file="seurat.all.markers.wilcox.res06.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.markers.wilcox.res06.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_alldata_filtered_monos.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_monos.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# read in markers to fix heatmaps 

#markers <- read.table(pipe('gsutil cat gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/seurat.all.markers.wilcox.res06.txt'))


top10 <- seurat.all.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


top10

top10$gene

#head(seurat.all.data.markers)

#seurat.all.data.markers %>% group_by(cluster) %>% top_n(10, p_val_adj) -> top10

DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_res06.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))
#ggsave("Exp_Heatmap_res02_BluYR.pdf", useDingbats = F)
#system(paste0("gsutil cp Exp_Heatmap_res02_BluYR.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


seurat.all.data.averages <- AverageExpression(seurat.all.data, return.seurat = TRUE)

DoHeatmap(seurat.all.data.averages,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_Avg_res06.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_Avg_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

DoHeatmap(seurat.all.data.averages,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))
ggsave("Exp_Heatmap_Avg_res06_BluYR.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_Avg_res06_BluYR.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


seurat.all.data <- BuildClusterTree(seurat.all.data, verbose = FALSE, reorder = TRUE)
tree <- PlotClusterTree(seurat.all.data)
ggsave("ClusterTree_res06.pdf", useDingbats = F)
system(paste0("gsutil cp ClusterTree_res02.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


head(seurat.all.data@meta.data)

seurat.all.factors <- seurat.all.data@meta.data$ArrayName
plate.table = table(seurat.all.factors, as.vector(seurat.all.data@meta.data$RNA_snn_res.0.6))
head(plate.table)

write.table(plate.table,file="cellspercluster_res06.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp cellspercluster_res06.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



FeaturePlot(seurat.all.data, features = c("CCR2","CD163","MSR1","CX3CR1","HLA-DR","CD14"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_Monocyte_FeaturePlot.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Monocyte_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



# 15 colors
colors.cluster.new <- c("#E86767","#BA2019","#BF6800","#FFF766","#B4C65C",
  "#029641","#48755C","#3FE0D8","#6A93CE","#767ADC",
  "#5C4A9D","#D9AAE5","#C872CF","#A44190","#C25183")

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("CCR2","CD163","CD68","MNDA","CX3CR1","HLA-DRA","CD14","CD1C","CST3")), pt.size = 0, col = colors.cluster.new)
ggsave("Monocyte_Violin_Plot_res06.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Monocyte_Violin_Plot_res06.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c(features = c("CTSS","IL1B","IL1RN","MT-CO1","CXCL5","FN1",
                                                   "FCER1A","RPS29","MALAT1","FCGR3A","IFI44L","FGFBP2","HBB","TUBB1","FCGR3B")),
                                        pt.size = 0, col = colors.cluster.new)
ggsave("Monocyte_Violin_Plot_res06_markers.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Monocyte_Violin_Plot_res06_markers.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


unique(colnames(seurat.all.data@meta.data))

## stacked barplot of clusters by donor, compartment, platform, and hours 
# donor (replicate)
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Replicate)) + geom_bar(position = "fill")
ggsave("StackedBarClusters_res06_replicate.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# compartment
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Compartment)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.compartment)
ggsave("StackedBarClusters_res06_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# platform
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Platform)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.platform)
ggsave("StackedBarClusters_res06_platform.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# hours
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Hours)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.hours)
ggsave("StackedBarClusters_res06_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# array 
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Hours)) + geom_bar(position = "fill") + scale_fill_manual(values=cols.hours)
ggsave("StackedBarClusters_res06_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_res06_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

# proper colors over array colored by consensus ID
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("StackedBarClusters_CellID_Array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_CellID_Array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)




# proper colors over array colored by cluster number
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new)
ggsave("StackedBarClusters_ClusterID_Array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "Compartment"
suerat.all.data.subset.BLD <- subset(seurat.all.data, idents = c("BLD"))
suerat.all.data.subset.HEF <- subset(seurat.all.data, idents = c("HEF"))


#blood doesn't have clusters: 4
colors.cluster.new.BLD  <- c("#E86767","#BA2019","#BF6800","#FFF766",
  "#029641","#48755C","#3FE0D8","#6A93CE","#767ADC",
  "#5C4A9D","#D9AAE5","#C872CF","#A44190","#C25183")


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.BLD) 
ggsave("StackedBarClusters_ClusterID_Array_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.BLD)
ggsave("StackedBarClusters_ClusterID_Hours_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# HEF doesnt have clusters: 7,10,12,13
colors.cluster.new.HEF <- c("#E86767","#BA2019","#BF6800","#FFF766","#B4C65C",
  "#029641","#48755C","#6A93CE","#767ADC",
  "#D9AAE5","#C25183")

# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.HEF)
ggsave("StackedBarClusters_ClusterID_Hours_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# proper colors over array colored by cluster number
ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=ArrayName, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.new.HEF)
ggsave("StackedBarClusters_ClusterID_Array_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Array_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


unique(seurat.all.data@meta.data$ArrayName)

Idents(object=seurat.all.data) <- "ArrayName"
suerat.all.data.subset.BLD.CTL <- subset(seurat.all.data, idents = c("DayF-500-BLD-00-Nova","DayF-500-BLD-01-Nova",
                                                                "DayF-500-BLD-02-Nova","DayF-500-BLD-03-Nova"))
suerat.all.data.subset.BLD.PT <- subset(seurat.all.data, idents = c("DayF-500-BLD-PT-Nova"))

# 15 colors
colors.cluster.CTL <- c("#E86767","#BA2019","#48755C","#3FE0D8",
                        "#6A93CE","#767ADC",
  "#5C4A9D","#D9AAE5","#C872CF","#A44190")

ggplot(suerat.all.data.subset.BLD.CTL@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.CTL)
ggsave("StackedBarClusters_ClusterID_Hours_BLD_Control.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_BLD_Control.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# 15 colors
colors.cluster.PT <- c("#E86767","#BA2019","#48755C","#3FE0D8","#6A93CE","#767ADC",
  "#5C4A9D","#D9AAE5","#A44190")

ggplot(suerat.all.data.subset.BLD.PT@meta.data, aes(x=Hours, fill=RNA_snn_res.0.6)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=colors.cluster.PT)
ggsave("StackedBarClusters_ClusterID_Hours_BLD_PT.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_ClusterID_Hours_BLD_PT.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# dotplot of receptor and ligands 

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
DotPlot(seurat.all.data, features = c("CCR9","CCR7","CXCR5","CXCR4","CCR10","CX3CR1","CCR8",
                                     "CCR5","CCR3","CCR2","CCR1","CXCR6","CXCR3","CCL1","CCL2","CCL18","MAP2K1",
                                     "CCL4","CCL7","CCL5","CXCL16","CXCL10")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("DotPlot.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp DotPlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)




gene.sets <- read.csv(pipe('gsutil cat gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/GeneSets_Mphage_Monos_DCs_MasterList.csv'))


head(gene.sets)
unique(colnames(gene.sets))

XUE_1 <- list(gene.sets$XUE_1)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_1,name = "XUE_1_", seed.use = 1988)


XUE_2 <- list(gene.sets$XUE_2)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_2,name = "XUE_2_", seed.use = 1988)

XUE_3 <- list(gene.sets$XUE_3)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_3,name = "XUE_3_", seed.use = 1988)

XUE_4 <- list(gene.sets$XUE_4)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_4,name = "XUE_4_", seed.use = 1988)

XUE_5 <- list(gene.sets$XUE_5)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_5,name = "XUE_5_", seed.use = 1988)

XUE_6 <- list(gene.sets$XUE_6)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_6,name = "XUE_6_", seed.use = 1988)

XUE_7 <- list(gene.sets$XUE_7)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_7,name = "XUE_7_", seed.use = 1988)

XUE_8 <- list(gene.sets$XUE_8)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_8,name = "XUE_8_", seed.use = 1988)

XUE_9 <- list(gene.sets$XUE_9)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_9,name = "XUE_9_", seed.use = 1988)

XUE_10 <- list(gene.sets$XUE_10)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_10,name = "XUE_10_", seed.use = 1988)

XUE_11 <- list(gene.sets$XUE_11)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_11,name = "XUE_11_", seed.use = 1988)

XUE_12 <- list(gene.sets$XUE_12)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_12,name = "XUE_12_", seed.use = 1988)

XUE_13 <- list(gene.sets$XUE_13)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_13,name = "XUE_13_", seed.use = 1988)

XUE_14 <- list(gene.sets$XUE_14)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_14,name = "XUE_14_", seed.use = 1988)

XUE_15 <- list(gene.sets$XUE_15)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_15,name = "XUE_15_", seed.use = 1988)

XUE_16 <- list(gene.sets$XUE_16)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_16,name = "XUE_16_", seed.use = 1988)

XUE_17 <- list(gene.sets$XUE_17)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_17,name = "XUE_17_", seed.use = 1988)

XUE_18 <- list(gene.sets$XUE_18)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_18,name = "XUE_18_", seed.use = 1988)

XUE_19 <- list(gene.sets$XUE_19)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_19,name = "XUE_19_", seed.use = 1988)

XUE_20 <- list(gene.sets$XUE_20)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_20,name = "XUE_20_", seed.use = 1988)

XUE_21 <- list(gene.sets$XUE_21)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_21,name = "XUE_21_", seed.use = 1988)

XUE_22 <- list(gene.sets$XUE_22)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_22,name = "XUE_22_", seed.use = 1988)

XUE_23 <- list(gene.sets$XUE_23)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_23,name = "XUE_23_", seed.use = 1988)

XUE_24 <- list(gene.sets$XUE_24)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_24,name = "XUE_24_", seed.use = 1988)

XUE_25 <- list(gene.sets$XUE_25)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_25,name = "XUE_25_", seed.use = 1988)

XUE_26 <- list(gene.sets$XUE_26)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_26,name = "XUE_26_", seed.use = 1988)

XUE_27 <- list(gene.sets$XUE_27)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_27,name = "XUE_27_", seed.use = 1988)

XUE_28 <- list(gene.sets$XUE_28)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_28,name = "XUE_28_", seed.use = 1988)

XUE_29 <- list(gene.sets$XUE_29)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_29,name = "XUE_29_", seed.use = 1988)

XUE_30 <- list(gene.sets$XUE_30)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_30,name = "XUE_30_", seed.use = 1988)

XUE_31 <- list(gene.sets$XUE_31)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_31,name = "XUE_31_", seed.use = 1988)

XUE_32 <- list(gene.sets$XUE_32)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_32,name = "XUE_32_", seed.use = 1988)

XUE_33 <- list(gene.sets$XUE_33)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_33,name = "XUE_33_", seed.use = 1988)

XUE_34 <- list(gene.sets$XUE_34)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_34,name = "XUE_34_", seed.use = 1988)

XUE_35 <- list(gene.sets$XUE_35)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_35,name = "XUE_35_", seed.use = 1988)

XUE_36 <- list(gene.sets$XUE_36)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_36,name = "XUE_36_", seed.use = 1988)

XUE_37 <- list(gene.sets$XUE_37)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_37,name = "XUE_37_", seed.use = 1988)

XUE_38 <- list(gene.sets$XUE_38)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_38,name = "XUE_38_", seed.use = 1988)

XUE_39 <- list(gene.sets$XUE_39)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_39,name = "XUE_39_", seed.use = 1988)

XUE_40 <- list(gene.sets$XUE_40)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_40,name = "XUE_40_", seed.use = 1988)

XUE_41 <- list(gene.sets$XUE_41)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_41,name = "XUE_41_", seed.use = 1988)

XUE_42 <- list(gene.sets$XUE_42)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_42,name = "XUE_42_", seed.use = 1988)

XUE_43 <- list(gene.sets$XUE_43)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_43,name = "XUE_43_", seed.use = 1988)

XUE_44 <- list(gene.sets$XUE_44)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_44,name = "XUE_44_", seed.use = 1988)

XUE_45 <- list(gene.sets$XUE_45)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_45,name = "XUE_45_", seed.use = 1988)

XUE_46 <- list(gene.sets$XUE_46)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_46,name = "XUE_46_", seed.use = 1988)

XUE_47 <- list(gene.sets$XUE_47)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_47,name = "XUE_47_", seed.use = 1988)

XUE_48 <- list(gene.sets$XUE_48)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_48,name = "XUE_48_", seed.use = 1988)

XUE_49 <- list(gene.sets$XUE_49)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = XUE_49,name = "XUE_49_", seed.use = 1988)

Monocyte_human <-list(gene.sets$Monocyte_human)
MyeloidDC_human <-list(gene.sets$MyeloidDC_human)
PlasmacytoidDC_human <-list(gene.sets$PlasmacytoidDC_human)

seurat.all.data <- AddModuleScore(object = seurat.all.data, features = Monocyte_human,name = "Monocyte_human", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = MyeloidDC_human,name = "MyeloidDC_human", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = PlasmacytoidDC_human,name = "PlasmacytoidDC_human", seed.use = 1988)

CD141_CLEC9A_Villani <-list(gene.sets$CD141_CLEC9A_Villani)
CD1C_A_Villani <-list(gene.sets$CD1C_A_Villani)
CD1C_B_Villani <-list(gene.sets$CD1C_B_Villani)
CD1Cminus_CD141minus_Villani <-list(gene.sets$CD1Cminus_CD141minus_Villani)
New_pop_Villani <-list(gene.sets$New_pop_Villani)
pDC_Villani <-list(gene.sets$pDC_Villani)

seurat.all.data <- AddModuleScore(object = seurat.all.data, features = CD141_CLEC9A_Villani,name = "CD141_CLEC9A_Villani", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = CD1C_A_Villani,name = "CD1C_A_Villani", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = CD1C_B_Villani,name = "CD1C_B_Villani", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = CD1Cminus_CD141minus_Villani,name = "CD1Cminus_CD141minus_Villani", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = New_pop_Villani,name = "New_pop_Villani", seed.use = 1988)
seurat.all.data <- AddModuleScore(object = seurat.all.data, features = pDC_Villani,name = "pDC_Villani", seed.use = 1988)


head(seurat.all.data@meta.data)

dim(seurat.all.data@meta.data)

FeaturePlot(seurat.all.data, features = c("XUE_10_1", "XUE_11_1","XUE_30_1","CD1C_B_Villani1"), pt.size = 0.5, reduction = "tsne", cols = c("gray100","darkslateblue"))
ggsave("tsne_Monocyte_Modules.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp tsne_Monocyte_Modules.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



unique(colnames(seurat.all.data@meta.data))

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c("XUE_10_1", "XUE_11_1","XUE_30_1","CD1C_B_Villani1"), pt.size = 0, y.max = 1.5, cols = colors.cluster.new)
ggsave("Monocyte_Violin_Plot_res06_GeneModules.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Monocyte_Violin_Plot_res06_GeneModules.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


seurat.all.factors <- seurat.all.data@meta.data$ArrayName
plate.table = table(seurat.all.factors, as.vector(seurat.all.data@meta.data$RNA_snn_res.0.6))
head(plate.table)

write.table(plate.table,file="cellspercluster_res06.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp cellspercluster_res06.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


metadata.forPCA <- seurat.all.data@meta.data
head(metadata.forPCA)

unique(colnames(metadata.forPCA))
dim(metadata.forPCA)

metadata.forPCA.annotation <- metadata.forPCA[,1:23]

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


module.scores <- metadata.forPCA[,24:81]

head(module.scores)

head(metadata.forPCA)
dim(metadata.forPCA)

row.annos <- cbind(metadata.forPCA.annotation$Hours,metadata.forPCA.annotation$Compartment)


row.annos <- cbind(row.annos, metadata.forPCA.annotation$Replicate)


head(row.annos)

row.annos <- cbind(row.annos,metadata.forPCA.annotation$RNA_snn_res.0.6)

row.annos <- as.data.frame(cbind(row.names(metadata.forPCA.annotation),row.annos))

row.names(row.annos) <- row.annos$V1
row.annos$V1 <- NULL
colnames(row.annos) <- c("Hours","Compartment","Replicate","Cluster")
head(row.annos)


unique(row.annos$Cluster) # 15 colors 
unique(row.annos$Hours) # 7 colors 
unique(row.annos$Compartment) # 2 colors 

row.annos$Hours <- paste("Hours", row.annos$Hours, sep="_")

row.annos$Cluster <- paste("Cluster", row.annos$Cluster, sep="_")

row.annos$Replicate <- paste("Replicate", row.annos$Replicate, sep="_")


head(row.annos)

row.annos[order(row.annos$Compartment),]
row.annos[order(row.annos$Hours),]

stim.annos <- read.csv(pipe('gsutil cat gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/module_stim_annotations.csv'))


head(stim.annos)

row.names(stim.annos) <- stim.annos$module
stim.annos$module <- NULL

right_ha = HeatmapAnnotation(df = stim.annos)


# 15 colors
colors.cluster.new <- c("#E86767","#BA2019","#BF6800","#FFF766","#B4C65C",
  "#029641",
                        "#48755C","#3FE0D8","#6A93CE","#767ADC","#5C4A9D",
                        "#D9AAE5","#C872CF","#A44190","#C25183")

colors_complex = list(Hours = c("Hours_050"="lightsteelblue1","Hours_066"="darkslategray1","Hours_073"="cadetblue2",
                                "Hours_089"="dodgerblue2","Hours_112"="blue1","Hours_137"="darkblue","Hours_500"="gold"),
                      Compartment = c("BLD"="brown","HEF"="darkgray"),
                      Cluster = c("Cluster_0"="#E86767","Cluster_1"="#BA2019","Cluster_2"="#BF6800",
                         "Cluster_3"="#FFF766","Cluster_4"="#B4C65C","Cluster_5"="#029641",
                         "Cluster_6"="#48755C","Cluster_7"="#3FE0D8","Cluster_8"="#6A93CE",
                         "Cluster_9"="#767ADC","Cluster_10"="#5C4A9D","Cluster_11"="#D9AAE5",
                         "Cluster_12"="#C872CF","Cluster_13"="#A44190","Cluster_14"="#C25183"),
                     Replicate = c("Replicate_R1"= "orange", "Replicate_R2" = "orange", "Replicate_00" = "gold", "Replicate_01" = "gold",
                                   "Replicate_02" = "gold", "Replicate_03" = "gold", "Replicate_PT" = "orange"))


top_ha = HeatmapAnnotation(df = row.annos, col = colors_complex)


Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE) + HeatmapAnnotation(df = stim.annos, which = "row", show_annotation_name = TRUE)


HA <- Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE) + HeatmapAnnotation(df = stim.annos, which = "row", show_annotation_name = TRUE)
rowsHA <- row_order(HA)
write.table(rowsHA,file="rowsHA.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp rowsHA.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

#save(HA,file="Monocyte_ModuleHeatmap.Robj")
#pdf("Monocyte_ModuleHeatmap.pdf")
#print(HA)
#dev.off()
#system(paste0("gsutil cp Monocyte_ModuleHeatmap.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


colnames(module.scores)

Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE) + HeatmapAnnotation(df = stim.annos, which = "row", show_annotation_name = TRUE)


help(row_order)

head(rowsHA)

head(stim.annos)

#HA <- Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE) + HeatmapAnnotation(df = stim.annos, which = "row")
#ggsave("Monocyte_ModuleHeatmap.pdf", useDingbats = F, height = 12, width = 16, units = "in")
#system(paste0("gsutil cp Monocyte_ModuleHeatmap.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE)

HA_NOL <- Heatmap(as.matrix(t(module.scores)), show_column_names = FALSE, top_annotation = top_ha,cluster_columns = FALSE, show_row_names = TRUE)
save(HA_NOL,file="Monocyte_ModuleHeatmap_norowannos.Robj")
pdf("Monocyte_ModuleHeatmap_norowannos.pdf")
print(HA_NOL)
dev.off()
system(paste0("gsutil cp Monocyte_ModuleHeatmap_norowannos.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_alldata_filtered_monos.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_monos.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


meta.data.with.modules <- seurat.all.data@meta.data
write.table(meta.data.with.modules,file="metadatawithmodules.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp metadatawithmodules.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


all.pca <- prcomp(module.scores)


ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, var.axes = FALSE,varname.adjust = 2) + scale_color_manual(values=colors.cluster.new)
ggsave("MonocyteBiplot_nolabels.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp MonocyteBiplot_nolabels.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2) + scale_color_manual(values=colors.cluster.new)
ggsave("MonocyteBiplot_labels.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp MonocyteBiplot_labels.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2,var.axes = FALSE, choices = 2:3) + scale_color_manual(values=colors.cluster.new)
ggsave("MonocyteBiplot_PCs2and3.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp MonocyteBiplot_PCs2and3.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


ggbiplot(all.pca, obs.scale = 1, var.scale = 1, groups = as.factor(row.annos$Cluster),
         circle = TRUE, varname.adjust = 2,var.axes = FALSE, choices = 3:4) + scale_color_manual(values=colors.cluster.new)
ggsave("MonocyteBiplot_PCs3and4.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp MonocyteBiplot_PCs3and4.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


var.loadings <- all.pca$rotation


write.table(var.loadings,file="variable_loadings.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp variable_loadings.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


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

head(var.loadings)

## make violin plots of modules: 

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("CD1C_B_Villani1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_CD1C_B_Villani1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_CD1C_B_Villani1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)




Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("CD1C_B_Villani1")), pt.size = 0, y.max = 1.2, cols = colors.cluster.new)


# 15 colors
#colors.cluster.new <- c(0"#E86767",1"#BA2019",2"#BF6800",3"#FFF766",4"#B4C65C",
#  5"#029641",6"#48755C",7"#3FE0D8",8"#6A93CE",9"#767ADC",
#  10"#5C4A9D",11"#D9AAE5",12"#C872CF",13"#A44190",14"#C25183")


colors.cluster.sort <- c("#BF6800","#BA2019","#E86767","#C872CF","#5C4A9D",
                        "#A44190","#FFF766","#3FE0D8","#D9AAE5","#C25183",
                        "#6A93CE","#48755C","#029641","#B4C65C","#767ADC")

Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
VlnPlot(seurat.all.data, features = c(features = c("CD1C_B_Villani1")), pt.size = 0, y.max = 1.2, sort = "increasing", cols = colors.cluster.sort)
ggsave("Tcell_Violin_Plot_res06_CD1C_B_Villani1_sorted.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_CD1C_B_Villani1_sorted.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


help(VlnPlot)

# PC1
VlnPlot(seurat.all.data, features = c(features = c("XUE_10_1")), pt.size = 0, y.max = 1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_XUE_10_1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_XUE_10_1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c(features = c("XUE_30_1")), pt.size = 0, y.max = 1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_XUE_30_1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_XUE_30_1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#PC2
VlnPlot(seurat.all.data, features = c(features = c("XUE_11_1")), pt.size = 0, y.max = 1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_XUE_11_1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_XUE_11_1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



VlnPlot(seurat.all.data, features = c(features = c("MyeloidDC_human1")), pt.size = 0, y.max = 1, cols = colors.cluster.new)
ggsave("Tcell_Violin_Plot_res06_MyeloidDC_human1.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Tcell_Violin_Plot_res06_MyeloidDC_human1.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pcs.matrix = all.pca$x[,c(1:2)]

help(save_plot)

pdf("AllModule_PCA.pdf")
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,1.5), ylim=c(-1,1), pch=19, cex=.3)
dev.off()
system(paste0("gsutil cp AllModule_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("CD1CBVillani_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,1.5), ylim=c(-1,1), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
dev.off()
#ggsave("CD1CBVillani_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp CD1CBVillani_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("XUE10_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,1.5), ylim=c(-1,1), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) # XUE 10 
dev.off()
#ggsave("XUE10_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp XUE10_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("XUE30_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,1.5), ylim=c(-1,1), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) # XUE 30
dev.off()
#ggsave("XUE10_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp XUE30_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("XUE11_PCA.pdf",useDingbats = F)
plot(pcs.matrix[,1], pcs.matrix[,2], xlim=c(-1,1.5), ylim=c(-1,1), pch=19, cex=.3, col = "cornsilk2")
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) # XUE 11
dev.off()
#ggsave("XUE11_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp XUE11_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


dim(var.loadings)

cells.use.c14 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_14")]
pdf("AllmoduleArrows_PCA.pdf",useDingbats = F)
plot(pcs.matrix[cells.use.c14,1], pcs.matrix[cells.use.c14,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 14, all mods", col="#e95dad")
arrows(0,0,var.loadings[1,][1],var.loadings[1,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[2,][1],var.loadings[2,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[3,][1],var.loadings[3,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[4,][1],var.loadings[4,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[5,][1],var.loadings[5,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[6,][1],var.loadings[6,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[7,][1],var.loadings[7,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[8,][1],var.loadings[8,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[9,][1],var.loadings[9,][2],col = "grey", code = 0)

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
arrows(0,0,var.loadings[29,][1],var.loadings[29,][2],col = "grey", code = 0)

arrows(0,0,var.loadings[31,][1],var.loadings[31,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[32,][1],var.loadings[32,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[33,][1],var.loadings[33,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[34,][1],var.loadings[34,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[35,][1],var.loadings[35,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[36,][1],var.loadings[36,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[37,][1],var.loadings[37,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[38,][1],var.loadings[38,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[39,][1],var.loadings[39,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[40,][1],var.loadings[30,][2],col = "grey", code = 0)

arrows(0,0,var.loadings[41,][1],var.loadings[41,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[42,][1],var.loadings[42,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[43,][1],var.loadings[43,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[44,][1],var.loadings[44,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[45,][1],var.loadings[45,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[46,][1],var.loadings[46,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[47,][1],var.loadings[47,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[48,][1],var.loadings[48,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[49,][1],var.loadings[49,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[50,][1],var.loadings[50,][2],col = "grey", code = 0)

arrows(0,0,var.loadings[51,][1],var.loadings[31,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[52,][1],var.loadings[32,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[53,][1],var.loadings[33,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[54,][1],var.loadings[34,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[56,][1],var.loadings[36,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[57,][1],var.loadings[37,][2],col = "grey", code = 0)
arrows(0,0,var.loadings[58,][1],var.loadings[38,][2],col = "grey", code = 0)

arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("AllmoduleArrows_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp AllmoduleArrows_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


cells.use.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF")]
pdf("Hematoma_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.hematoma,1], pcs.matrix[cells.use.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="hematoma", col="darkgray")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Hematoma_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Hematoma_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_PCA.pdf", useDingbats = F)
cells.use.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD")]
plot(pcs.matrix[cells.use.blood,1], pcs.matrix[cells.use.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="blood", col="brown")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Patient_PCA.pdf", useDingbats = F)
cells.use.PT = rownames(row.annos)[which(row.annos$Replicate=="Replicate_R1" | row.annos$Replicate=="Replicate_R2" | row.annos$Replicate=="Replicate_PT")]
plot(pcs.matrix[cells.use.PT,1], pcs.matrix[cells.use.PT,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="Patient", col="orange")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
system(paste0("gsutil cp Patient_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Control_PCA.pdf", useDingbats = F)
cells.use.controls = rownames(row.annos)[which(row.annos$Replicate=="Replicate_00" | row.annos$Replicate=="Replicate_01" | row.annos$Replicate=="Replicate_02" | row.annos$Replicate=="Replicate_03")]
plot(pcs.matrix[cells.use.controls,1], pcs.matrix[cells.use.controls,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="Controls", col="gold")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
system(paste0("gsutil cp Control_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


head(row.annos)
unique(row.annos$Hours)

# cluster 0
pdf("Cluster0_PCA.pdf", useDingbats = F)
cells.use.c0 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_0")]
plot(pcs.matrix[cells.use.c0,1], pcs.matrix[cells.use.c0,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 0", col="#E86767")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster0_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster0_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 1
cells.use.c1 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_1")]
pdf("Cluster1_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c1,1], pcs.matrix[cells.use.c1,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 1", col="#BA2019")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster1_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster1_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 2
cells.use.c2 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_2")]
pdf("Cluster2_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c2,1], pcs.matrix[cells.use.c2,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 2", col="#BF6800")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster2_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster2_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 3
cells.use.c3 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_3")]
pdf("Cluster3_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c3,1], pcs.matrix[cells.use.c3,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 3", col="#FFF766")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster3_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster3_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 4
cells.use.c4 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_4")]
pdf("Cluster4_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c4,1], pcs.matrix[cells.use.c4,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 4", col="#B4C65C")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster4_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster4_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 5
cells.use.c5 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_5")]
pdf("Cluster5_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c5,1], pcs.matrix[cells.use.c5,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 5", col="#029641")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster5_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster5_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 6
cells.use.c6 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_6")]
pdf("Cluster6_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c6,1], pcs.matrix[cells.use.c6,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 6", col="#48755C")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster6_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster6_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 7
cells.use.c7 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_7")]
pdf("Cluster7_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c7,1], pcs.matrix[cells.use.c7,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 7", col="#3FE0D8")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster7_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster7_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 8
cells.use.c8 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_8")]
pdf("Cluster8_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c8,1], pcs.matrix[cells.use.c8,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 8", col="#6A93CE")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster8_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster8_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 9
cells.use.c9 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_9")]
pdf("Cluster9_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c9,1], pcs.matrix[cells.use.c9,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 9", col="#767ADC")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster9_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster9_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 10
cells.use.10 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_10")]
pdf("Cluster10_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.10,1], pcs.matrix[cells.use.10,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 10", col="#5C4A9D")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster10_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster10_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 11
cells.use.c11 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_11")]
pdf("Cluster11_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c11,1], pcs.matrix[cells.use.c11,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 11", col="#D9AAE5")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster11_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster11_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#Cluster = c("Cluster_0"="#E86767","Cluster_1"="#BA2019","Cluster_2"="#BF6800",
#                         "Cluster_3"="#FFF766","Cluster_4"="#B4C65C","Cluster_5"="#029641",
#                         "Cluster_6"="#48755C","Cluster_7"="#3FE0D8","Cluster_8"="#6A93CE",
#                         "Cluster_9"="#767ADC","Cluster_10"="#5C4A9D","Cluster_11"="#D9AAE5",
#                         "Cluster_12"="#C872CF","Cluster_13"="#A44190","Cluster_14"="#C25183")

# cluster 12
cells.use.c12 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_12")]
pdf("Cluster12_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c12,1], pcs.matrix[cells.use.c12,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 12", col="#C872CF")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster12_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster12_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 13
cells.use.c13 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_13")]
pdf("Cluster13_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c13,1], pcs.matrix[cells.use.c13,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 13", col="#A44190")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster13_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster13_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


# cluster 14
cells.use.c14 = rownames(row.annos)[which(row.annos$Cluster=="Cluster_14")]
pdf("Cluster14_PCA.pdf", useDingbats = F)
plot(pcs.matrix[cells.use.c14,1], pcs.matrix[cells.use.c14,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="cluster 14", col="#C25183")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Cluster14_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Cluster14_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_50hrs_PCA.pdf", useDingbats = F)
cells.use.50hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_050")]
plot(pcs.matrix[cells.use.50hrs.blood,1], pcs.matrix[cells.use.50hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="50 hours, blood", col="lightsteelblue1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_50hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_50hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_66hrs_PCA.pdf", useDingbats = F)
cells.use.66hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_066")]
plot(pcs.matrix[cells.use.66hrs.blood,1], pcs.matrix[cells.use.66hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="66 hours, blood", col="darkslategray1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_66hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_66hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_73hrs_PCA.pdf", useDingbats = F)
cells.use.73hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_073")]
plot(pcs.matrix[cells.use.73hrs.blood,1], pcs.matrix[cells.use.73hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="73 hours, blood", col="cadetblue2")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_73hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_73hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_89hrs_PCA.pdf", useDingbats = F)
cells.use.89hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_089")]
plot(pcs.matrix[cells.use.89hrs.blood,1], pcs.matrix[cells.use.89hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="89 hours, blood", col="dodgerblue2")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_89hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_89hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_112hrs_PCA.pdf", useDingbats = F)
cells.use.112hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_112")]
plot(pcs.matrix[cells.use.112hrs.blood,1], pcs.matrix[cells.use.112hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="112 hours, blood", col="blue1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_112hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_112hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_137hrs_PCA.pdf", useDingbats = F)
cells.use.137hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_137")]
plot(pcs.matrix[cells.use.137hrs.blood,1], pcs.matrix[cells.use.137hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="137 hours, blood", col="darkblue")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_137hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_137hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("Blood_500hrs_PCA.pdf", useDingbats = F)
cells.use.500hrs.blood = rownames(row.annos)[which(row.annos$Compartment=="BLD" & row.annos$Hours=="Hours_500")]
plot(pcs.matrix[cells.use.500hrs.blood,1], pcs.matrix[cells.use.500hrs.blood,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="500 hours, blood", col="gold")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("Blood_500hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Blood_500hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


pdf("hematoma_50hrs_PCA.pdf", useDingbats = F)
cells.use.50hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_050")]
plot(pcs.matrix[cells.use.50hrs.hematoma,1], pcs.matrix[cells.use.50hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="50 hours, hematoma", col="lightsteelblue1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_50hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_50hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

pdf("hematoma_66hrs_PCA.pdf", useDingbats = F)
cells.use.66hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_066")]
plot(pcs.matrix[cells.use.66hrs.hematoma,1], pcs.matrix[cells.use.66hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="66 hours, hematoma", col="darkslategray1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_66hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_66hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

pdf("hematoma_73hrs_PCA.pdf", useDingbats = F)
cells.use.73hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_073")]
plot(pcs.matrix[cells.use.73hrs.hematoma,1], pcs.matrix[cells.use.73hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="73 hours, hematoma", col="cadetblue2")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_73hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_73hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

pdf("hematoma_89hrs_PCA.pdf", useDingbats = F)
cells.use.89hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_089")]
plot(pcs.matrix[cells.use.89hrs.hematoma,1], pcs.matrix[cells.use.89hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="89 hours, hematoma", col="dodgerblue2")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_89hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_89hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

pdf("hematoma_112hrs_PCA.pdf", useDingbats = F)
cells.use.112hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_112")]
plot(pcs.matrix[cells.use.112hrs.hematoma,1], pcs.matrix[cells.use.112hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="112 hours, hematoma", col="blue1")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_112hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_112hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

pdf("hematoma_137hrs_PCA.pdf", useDingbats = F)
cells.use.137hrs.hematoma = rownames(row.annos)[which(row.annos$Compartment=="HEF" & row.annos$Hours=="Hours_137")]
plot(pcs.matrix[cells.use.137hrs.hematoma,1], pcs.matrix[cells.use.137hrs.hematoma,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.7, main="137 hours, hematoma", col="darkblue")
arrows(0,0,var.loadings[55,][1],var.loadings[55,][2]) # CD1C_B Villani 
arrows(0,0,var.loadings[10,][1],var.loadings[10,][2]) #XUE 10 
arrows(0,0,var.loadings[30,][1],var.loadings[30,][2]) #XUE 30
arrows(0,0,var.loadings[11,][1],var.loadings[11,][2]) #XUE 11
dev.off()
#ggsave("hematoma_137hrs_PCA.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp hematoma_137hrs_PCA.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



head(seurat.all.data@meta.data)
colnames(seurat.all.data@meta.data)

# make a new column in Seurat object and rename clusters there before extracting 
Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
seurat.all.data <- RenameIdents(seurat.all.data, 
`0` = "cluster0", `1` = "cluster1", `2` = "cluster2", 
`3` = "cluster3", `4` = "cluster4", `5` = "cluster5", 
`6` = "cluster6", `7` = "cluster7", `8` = "cluster8", 
`9` = "cluster9", `10` = "cluster10", `11` = "cluster11", 
`12` = "cluster12", `13` = "cluster13", `14` = "cluster14")

seurat.all.data <- StashIdent(object = seurat.all.data, save.name = 'myeloidsubclusters')

colnames(seurat.all.data@meta.data)

data.for.anova <- seurat.all.data@meta.data

colnames(data.for.anova)
dim(data.for.anova)

data.for.anova <- data.for.anova[,78:82]

head(data.for.anova)

data.for.anova$CD1Cminus_CD141minus_Villani1 <- NULL
data.for.anova$New_pop_Villani1 <- NULL
data.for.anova$pDC_Villani1 <- NULL


colnames(data.for.anova) <- c("CD1C_B_Villani1_Score", "Cluster")

#rownames(data.for.anova) <- NULL

#CD1Cscore_cleaned <- as.data.frame(cbind(seurat.all.data@meta.data$CD1C_B_Villani1,
#                                         seurat.all.data@meta.data$RNA_snn_res.0.6))

#CD1Cscore_cleaned <- as.data.frame(cbind(data.for.anova$CD1C_B_Villani1,
#                            data.for.anova$myeloidsubclusters))


#head(CD1Cscore_cleaned)

#colnames(CD1Cscore_cleaned) <- c("CD1C_B_Villani1_Score", "Cluster")

#head(CD1Cscore_cleaned)

#unique(CD1Cscore_cleaned$Cluster)

#CD1Cscore_cleaned$Cluster <- paste("Cluster", CD1Cscore_cleaned$Cluster, sep ="_")

head(data.for.anova)
class(data.for.anova)

class(data.for.anova)
dim(data.for.anova)

unique(data.for.anova$Cluster)

CD1Cscore_cleaned <- data.for.anova #reassign

hist(CD1Cscore_cleaned$CD1C_B_Villani1_Score)

plotmeans(CD1Cscore_cleaned$CD1C_B_Villani1_Score~CD1Cscore_cleaned$Cluster, digits=2, ccol="red",mean.labels=T, main="plot of means")
ggsave("meanplot_ANOVA_CD1C_BScores.pdf", useDingbats = F, height = 12, width = 12, units = "in")
system(paste0("gsutil cp meanplot_ANOVA_CD1C_BScores.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)



aov_CD1Cscore <- aov(CD1Cscore_cleaned$CD1C_B_Villani1_Score ~ CD1Cscore_cleaned$Cluster)
summary(aov_CD1Cscore)

# 1.1 Are variance homogeneous? 
# plot residuals - only three outliers ID'd. Overall, suggests homoegeous residuals, run Tukey post test for multiple comparisons. 
plot(aov_CD1Cscore, 1)


# 1.2 Are variance homogeneous? 
# use bartlett test to determine if variance is the same for each group. 
# Null = variances are the same. Alternative = variance not the same. 
# P value is very small, so we reject the null, and variance is not the same for each group. 
bartlett.test(CD1C_B_Villani1_Score ~ Cluster, data = CD1Cscore_cleaned)


# 2. check for normality - majority of points fall along reference line, but tails deviate, so not normal. 
plot(aov_CD1Cscore, 2)


# confirm with Shapiro-Wilk test (note: dataset too large to run with this test)
aov_residuals <- residuals(object = aov_CD1Cscore)
shapiro.test(x = aov_residuals)


tuk_CD1Cscore <- TukeyHSD(aov_CD1Cscore)
results_tuk <- as.data.frame(tuk_CD1Cscore$`CD1Cscore_cleaned$Cluster`)


tuk_CD1Cscore

plot(tuk_CD1Cscore)
ggsave("Famwiseconfidencelevel_ANOVA_CD1CScore.pdf", useDingbats = F, height = 12, width = 12, units = "in")
system(paste0("gsutil cp Famwiseconfidencelevel_ANOVA_CD1CScore.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


write.table(results_tuk,file="ANOVA_Tukey_CD1CBScores.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp ANOVA_Tukey_CD1CBScores.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


### Run non-parametric alternative with Kruskal-Wallis then use pairwise wilcox to fiture out which pairs of groups are different 


KW_CD1Cscore <- kruskal.test(CD1Cscore_cleaned$CD1C_B_Villani1_Score ~ CD1Cscore_cleaned$Cluster)

KW_CD1Cscore # there are significant differences between scores in clusters

### pairwise compairisons

pairwise_wilcox <- pairwise.wilcox.test(CD1Cscore_cleaned$CD1C_B_Villani1_Score, CD1Cscore_cleaned$Cluster,
                 p.adjust.method = "BH")

pairwise_wilcox

class(pairwise_wilcox)

res_ranksum <- pairwise_wilcox$p.value

write.table(res_ranksum,file="WilcoxRankSum_CD1CBScores.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp WilcoxRankSum_CD1CBScores.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


hist(pairwise_wilcox$p.value)

colnames(seurat.all.data@meta.data)
unique(seurat.all.data@meta.data$idents_res2_consensus)

Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.monosonly <- subset(seurat.all.data, idents = c("CD1C+ DCs"), invert = TRUE)

unique(suerat.all.data.monosonly@meta.data$idents_res2_consensus)

Idents(object=suerat.all.data.monosonly) <- "Hours"
suerat.all.data.monosonly.050 <- subset(suerat.all.data.monosonly, idents = c("050"))
suerat.all.data.monosonly.066 <- subset(suerat.all.data.monosonly, idents = c("066"))
suerat.all.data.monosonly.073 <- subset(suerat.all.data.monosonly, idents = c("073"))
suerat.all.data.monosonly.089 <- subset(suerat.all.data.monosonly, idents = c("089"))
suerat.all.data.monosonly.112 <- subset(suerat.all.data.monosonly, idents = c("112"))
suerat.all.data.monosonly.137 <- subset(suerat.all.data.monosonly, idents = c("137"))

Idents(object=suerat.all.data.monosonly.050) <- "Compartment"
suerat.all.data.monosonly.050.BLD <- subset(suerat.all.data.monosonly.050, idents = c("BLD"))
suerat.all.data.monosonly.050.HEF <- subset(suerat.all.data.monosonly.050, idents = c("HEF"))
Idents(object=suerat.all.data.monosonly.066) <- "Compartment"
suerat.all.data.monosonly.066.BLD <- subset(suerat.all.data.monosonly.066, idents = c("BLD"))
suerat.all.data.monosonly.066.HEF <- subset(suerat.all.data.monosonly.066, idents = c("HEF"))
Idents(object=suerat.all.data.monosonly.073) <- "Compartment"
suerat.all.data.monosonly.073.BLD <- subset(suerat.all.data.monosonly.073, idents = c("BLD"))
suerat.all.data.monosonly.073.HEF <- subset(suerat.all.data.monosonly.073, idents = c("HEF"))
Idents(object=suerat.all.data.monosonly.089) <- "Compartment"
suerat.all.data.monosonly.089.BLD <- subset(suerat.all.data.monosonly.089, idents = c("BLD"))
suerat.all.data.monosonly.089.HEF <- subset(suerat.all.data.monosonly.089, idents = c("HEF"))
Idents(object=suerat.all.data.monosonly.112) <- "Compartment"
suerat.all.data.monosonly.112.BLD <- subset(suerat.all.data.monosonly.112, idents = c("BLD"))
suerat.all.data.monosonly.112.HEF <- subset(suerat.all.data.monosonly.112, idents = c("HEF"))
Idents(object=suerat.all.data.monosonly.137) <- "Compartment"
suerat.all.data.monosonly.137.BLD <- subset(suerat.all.data.monosonly.137, idents = c("BLD"))
suerat.all.data.monosonly.137.HEF <- subset(suerat.all.data.monosonly.137, idents = c("HEF"))

# 50 
count.data.050.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.050.BLD, slot = "counts"))
count.data.050.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.050.HEF, slot = "counts"))

count.data.050.sum <- cbind(rowSums(count.data.050.BLD),rowSums(count.data.050.HEF))
colnames(count.data.050.sum) <- c("pseudobulk_monos_BLD_50hrs","pseudobulk_monos_HEF_50hrs")

write.table(count.data.050.sum,file="Pseudobulk_050hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_050hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)

head(count.data.050.sum)

# 66
count.data.066.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.066.BLD, slot = "counts"))
count.data.066.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.066.HEF, slot = "counts"))

count.data.066.sum <- cbind(rowSums(count.data.066.BLD),rowSums(count.data.066.HEF))
colnames(count.data.066.sum) <- c("pseudobulk_monos_BLD_66hrs","pseudobulk_monos_HEF_66hrs")

write.table(count.data.066.sum,file="Pseudobulk_066hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_066hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#73
count.data.073.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.073.BLD, slot = "counts"))
count.data.073.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.073.HEF, slot = "counts"))

count.data.073.sum <- cbind(rowSums(count.data.073.BLD),rowSums(count.data.073.HEF))
colnames(count.data.073.sum) <- c("pseudobulk_monos_BLD_73hrs","pseudobulk_monos_HEF_73hrs")

write.table(count.data.073.sum,file="Pseudobulk_073hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_073hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#89
count.data.089.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.089.BLD, slot = "counts"))
count.data.089.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.089.HEF, slot = "counts"))

count.data.089.sum <- cbind(rowSums(count.data.089.BLD),rowSums(count.data.089.HEF))
colnames(count.data.089.sum) <- c("pseudobulk_monos_BLD_89hrs","pseudobulk_monos_HEF_89hrs")

write.table(count.data.089.sum,file="Pseudobulk_089hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_089hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#112
count.data.112.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.112.BLD, slot = "counts"))
count.data.112.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.112.HEF, slot = "counts"))

count.data.112.sum <- cbind(rowSums(count.data.112.BLD),rowSums(count.data.112.HEF))
colnames(count.data.112.sum) <- c("pseudobulk_monos_BLD_112hrs","pseudobulk_monos_HEF_112hrs")

write.table(count.data.112.sum,file="Pseudobulk_112hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_112hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


#137
count.data.137.BLD <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.137.BLD, slot = "counts"))
count.data.137.HEF <- as.data.frame(GetAssayData(object = suerat.all.data.monosonly.137.HEF, slot = "counts"))

count.data.137.sum <- cbind(rowSums(count.data.137.BLD),rowSums(count.data.137.HEF))
colnames(count.data.137.sum) <- c("pseudobulk_monos_BLD_137hrs","pseudobulk_monos_HEF_137hrs")

write.table(count.data.137.sum,file="Pseudobulk_137hrs_monos.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp Pseudobulk_137hrs_monos.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster6 <- subset(seurat.all.data, idents = c("6"))


seurat.factors.cluster6 <- suerat.all.data.subset.cluster6@meta.data$RNA_snn_res.0.6
plate.table.cluster6 = table(seurat.factors.cluster6, as.vector(suerat.all.data.subset.cluster6@meta.data$Compartment))
plate.table.cluster6


Idents(object=suerat.all.data.subset.cluster6) <- "Compartment"

seurat.all.data.markers.cluster6 <- FindAllMarkers(suerat.all.data.subset.cluster6, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 87, random.seed = 1988)

write.table(seurat.all.data.markers.cluster6,file="seurat.all.data.markers.cluster6.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster6.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


top10.cluster6 <- seurat.all.data.markers.cluster6 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster6,features=top10.cluster6$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster6.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster6.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster9 <- subset(seurat.all.data, idents = c("9"))


seurat.factors.cluster9 <- suerat.all.data.subset.cluster9@meta.data$RNA_snn_res.0.6
plate.table.cluster9 = table(seurat.factors.cluster9, as.vector(suerat.all.data.subset.cluster9@meta.data$Compartment))
plate.table.cluster9


Idents(object=suerat.all.data.subset.cluster9) <- "Compartment"

seurat.all.data.markers.cluster9 <- FindAllMarkers(suerat.all.data.subset.cluster9, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 27, random.seed = 1988)

write.table(seurat.all.data.markers.cluster9,file="seurat.all.data.markers.cluster9.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster9.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


top10.cluster9 <- seurat.all.data.markers.cluster9 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster9,features=top10.cluster9$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster9.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster9.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


Idents(object=seurat.all.data) <- "RNA_snn_res.0.6"
suerat.all.data.subset.cluster14 <- subset(seurat.all.data, idents = c("14"))


seurat.factors.cluster14 <- suerat.all.data.subset.cluster14@meta.data$RNA_snn_res.0.6
plate.table.cluster14 = table(seurat.factors.cluster14, as.vector(suerat.all.data.subset.cluster14@meta.data$Compartment))
plate.table.cluster14


Idents(object=suerat.all.data.subset.cluster14) <- "Compartment"

seurat.all.data.markers.cluster14 <- FindAllMarkers(suerat.all.data.subset.cluster14, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5, max.cells.per.ident = 11, random.seed = 1988)

write.table(seurat.all.data.markers.cluster14,file="seurat.all.data.markers.cluster14.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.data.markers.cluster14.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)


top10.cluster14 <- seurat.all.data.markers.cluster14 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(suerat.all.data.subset.cluster14,features=top10.cluster14$gene, raster = FALSE) + theme(text = element_text(size = 8))
ggsave("Exp_Heatmap_cluster14.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_cluster14.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200508_MonocytesAnalysis/"), intern=TRUE)






