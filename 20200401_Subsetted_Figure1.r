library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library('data.table')
library(stringr)
library(dplyr)
library(RColorBrewer)

project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
project
workspace
bucket

list.files(path = ".")

#RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200329_MultipletsRemoved_Figure1/*.rds"
#import updated obj to continue analyses 

RDS_Files <- "gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/*.rds"

system(sprintf("gsutil ls %s", RDS_Files), intern=T)

system("mkdir $(pwd)/RDS_Files")
system(sprintf("gsutil -m cp %s $(pwd)/RDS_Files/", RDS_Files))

list.files(path = "./RDS_Files")

#seurat.all.data <- readRDS("./RDS_Files/seuratobj_subsetted_cells.rds") #original import 
seurat.all.data <- readRDS("./RDS_Files/seuratobj_alldata_filtered.rds") #import for re-analysis

head(seurat.all.data@meta.data)

unique(seurat.all.data@meta.data$ArrayName)

#remove original clustering 
#seurat.all.data$RNA_snn_res.0.5 = NULL
#seurat.all.data$RNA_snn_res.0.6 = NULL
#seurat.all.data$RNA_snn_res.0.8 = NULL
#seurat.all.data$RNA_snn_res.1 = NULL
#seurat.all.data$RNA_snn_res.1.5 = NULL
#seurat.all.data$firstpassres2 = seurat.all.data$RNA_snn_res.2
#seurat.all.data$RNA_snn_res.2 = NULL

dim(seurat.all.data@meta.data)

# colors 
cols.hours = c("lightsteelblue1","darkslategray1","cadetblue2","dodgerblue2","blue1","darkblue","gold")
cols.compartment = c("brown","darkgray")
cols.platform = c("rosybrown","peru")

Idents(object=seurat.all.data) <- "Hours"

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1, cols = cols.hours)
ggsave("nFeatureRNA_violin_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp nFeatureRNA_violin_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1,cols = cols.hours)
ggsave("nCountRNA_violin_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp nCountRNA_violin_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1,cols = cols.hours)
ggsave("ptMT_violin_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp ptMT_violin_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, cols = cols.hours)
ggsave("techviolin_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp techviolin_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


Idents(object=seurat.all.data) <- "Platform"

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1, cols = cols.platform)
ggsave("nFeatureRNA_violin_Platform.pdf", useDingbats = F)
system(paste0("gsutil cp nFeatureRNA_violin_Platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1, cols = cols.platform)
ggsave("nCountRNA_violin_Platform.pdf", useDingbats = F)
system(paste0("gsutil cp nCountRNA_violin_Platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1, cols = cols.platform)
ggsave("ptMT_violin_Platform.pdf", useDingbats = F)
system(paste0("gsutil cp ptMT_violin_Platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, cols = cols.platform)
ggsave("techviolin_Platform.pdf", useDingbats = F)
system(paste0("gsutil cp techviolin_Platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



Idents(object=seurat.all.data) <- "Compartment"

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1, cols = cols.compartment)
ggsave("nFeatureRNA_violin_Compartment.pdf", useDingbats = F)
system(paste0("gsutil cp nFeatureRNA_violin_Compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1, cols = cols.compartment)
ggsave("nCountRNA_violin_Compartment.pdf", useDingbats = F)
system(paste0("gsutil cp nCountRNA_violin_Compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1, cols = cols.compartment)
ggsave("ptMT_violin_Compartment.pdf", useDingbats = F)
system(paste0("gsutil cp ptMT_violin_Compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, cols = cols.compartment)
ggsave("techviolin_Compartment.pdf", useDingbats = F)
system(paste0("gsutil cp techviolin_Compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

Idents(object=seurat.all.data) <- "ArrayName"

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1)
ggsave("nFeatureRNA_violin_Replicate.pdf", useDingbats = F, height = 7, width = 10, units = "in")
system(paste0("gsutil cp nFeatureRNA_violin_Replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1)
ggsave("nCountRNA_violin_Replicate.pdf", useDingbats = F, height = 7, width = 10, units = "in")
system(paste0("gsutil cp nCountRNA_violin_Replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1)
ggsave("ptMT_violin_Replicate.pdf", useDingbats = F, height = 7, width = 10, units = "in")
system(paste0("gsutil cp ptMT_violin_Replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
ggsave("techviolin_Replicate.pdf", useDingbats = F, height = 7, width = 10, units = "in")
system(paste0("gsutil cp techviolin_Replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)





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
system(paste0("gsutil cp Variable_Genes_Merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# scale and PCA, can regress stuff out here too (do not regress for now)
seurat.all.data <- ScaleData(seurat.all.data, features = rownames(seurat.all.data))
seurat.all.data <- RunPCA(seurat.all.data, features = VariableFeatures(object = seurat.all.data))

ElbowPlot(seurat.all.data)
ggsave("Elbow_merge.pdf", useDingbats = F)
system(paste0("gsutil cp Elbow_merge.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


## Jackstraw sig PCs (takes a really long time to run - all 20 PCs are sig!)
#seurat.all.data <- JackStraw(seurat.all.data, num.replicate = 100)
#seurat.all.data <- ScoreJackStraw(seurat.all.data, dims = 1:20)
#JackStrawPlot(seurat.all.data, dims = 1:20)
ggsave("JackstrawPlot.pdf", useDingbats = F)
system(paste0("gsutil cp JackstrawPlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

################# UMAP and clustering
seurat.all.data <- RunUMAP(seurat.all.data, reduction = "pca", dims = 1:17)
seurat.all.data <- RunTSNE(seurat.all.data, reduction = "pca", dims = 1:20)
seurat.all.data <- FindNeighbors(seurat.all.data, reduction = "pca", dims = 1:20)

# QC feature plots
FeaturePlot(seurat.all.data, features = c("percent.mt"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_percentmito.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_percentmito.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

FeaturePlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 1, cols = c("gray100","firebrick1"), reduction = "tsne")
ggsave("tsne_RNAcount.pdf", useDingbats = F)
system(paste0("gsutil cp tsne_RNAcount.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# colored by other features
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform)
ggsave("UMAP_platform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform")
ggsave("UMAP_platform_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment)
ggsave("UMAP_compartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment")
ggsave("UMAP_compartment_split.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8)
ggsave("UMAP_arrayname.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours)
ggsave("UMAP_hours.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours")
ggsave("UMAP_hours_splitbyhrs.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate")
ggsave("UMAP_hours_splitbyreplicate.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform")
ggsave("UMAP_hours_splitbyplatform.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment")
ggsave("UMAP_hours_splitbycompartment.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# colored by other features (no legend)
# platform 
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform) + NoLegend()
ggsave("UMAP_platform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Platform", pt.size = 0.8, cols = cols.platform, split.by = "Platform") + NoLegend()
ggsave("UMAP_platform_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_platform_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#compartment
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment) + NoLegend()
ggsave("UMAP_compartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Compartment", pt.size = 0.8, cols = cols.compartment, split.by = "Compartment") + NoLegend()
ggsave("UMAP_compartment_split_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_compartment_split_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#Array
DimPlot(seurat.all.data, reduction = "tsne", group.by = "ArrayName", pt.size = 0.8) + NoLegend()
ggsave("UMAP_arrayname_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_arrayname_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# hours
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours) + NoLegend()
ggsave("UMAP_hours_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols =cols.hours, split.by = "Hours") + NoLegend()
ggsave("UMAP_hours_splitbyhrs_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Replicate") + NoLegend()
ggsave("UMAP_hours_splitbyreplicate_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Platform") + NoLegend()
ggsave("UMAP_hours_splitbyplatform_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyreplicate_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Compartment") + NoLegend()
ggsave("UMAP_hours_splitbycompartment_NoL.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbycompartment_NoL.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


## DC/Mono/Mac feature plots
FeaturePlot(seurat.all.data, features = c("CD68", "MSR1", "CD14","ITGAM", 
                                          "CD1C", "MNDA", "CST3","HLA-DRA"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_DCmonoMac_FeaturePlot.pdf", useDingbats = F, height = 9, width = 9, units = "in")
system(paste0("gsutil cp tsne_DCmonoMac_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

## T cell feature plots 
FeaturePlot(seurat.all.data, features = c("TRAC", "CD4", "CD8A", "IL7R",
                                          "TRA","TRBC2","CD3E","CD3G"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_Tcell_FeaturePlot.pdf", useDingbats = F, height = 9, width = 9, units = "in")
system(paste0("gsutil cp tsne_Tcell_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


## B cells and NK cells
FeaturePlot(seurat.all.data, features = c("GNLY","NKG7","GZMA","MS4A7", "CD19", "MS4A1"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_BcellNKcell_FeaturePlot.pdf", useDingbats = F, height = 9, width = 9, units = "in")
system(paste0("gsutil cp tsne_BcellNKcell_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


FeaturePlot(seurat.all.data, features = c("KIF5A", "HBB", "HBA2"), pt.size = 0.5, reduction = "tsne")
ggsave("tsne_NeuronandRBC_FeaturePlot.pdf", useDingbats = F, height = 9, width = 9, units = "in")
system(paste0("gsutil cp tsne_NeuronandRBC_FeaturePlot.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# cluster IDs
seurat.all.data <- FindClusters(seurat.all.data, resolution = 2)


# all together
DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5)  
ggsave("tsne_clusterID.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("tsne_clusterID_no_legend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterIDno_legend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#split by compartment
DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment")  
ggsave("tsne_clusterID_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Compartment") + NoLegend()
ggsave("tsne_clusterID_no_legend_compartmentsplit.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_no_legend_compartmentsplit.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)




# Violin plots of key lineage markers for current cluster IDs (not run for 0.3)
Idents(object=seurat.all.data) <- "seurat_clusters"
VlnPlot(seurat.all.data, features = c("GNLY","NKG7","GZMA","CD74","CD19","IGHM"), pt.size = 0)
ggsave("NKandBcell_Violin_ClusterID.pdf", useDingbats = F, height = 8, width = 12, units = "in")
system(paste0("gsutil cp NKandBcell_Violin_ClusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("CD68", "MSR1", "FCGR3A", "CD14", "ITGAM",
                                      "CD1C", "CST3", "MS4A7","MNDA"), pt.size = 0)
ggsave("DCMonoMac_Violin_ClusterID.pdf", useDingbats = F, height = 8, width = 12, units = "in")
system(paste0("gsutil cp DCMonoMac_Violin_ClusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("PTPRC", "TRAC", "CD4", "CD8A", "IL7R",
                                       "TRDC","TRBC2","CD3E","CD3G"), pt.size = 0)
ggsave("Tcell_Violin_ClusterID.pdf", useDingbats = F, height = 8, width = 12, units = "in")
system(paste0("gsutil cp Tcell_Violin_ClusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("KIF5A", "HBB", "HBA2"), pt.size = 0)
ggsave("NeuronandRBC_Violin_ClusterID.pdf", useDingbats = F, height = 8, width = 12, units = "in")
system(paste0("gsutil cp NeuronandRBC_Violin_ClusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



seurat.all.data.markers <- FindAllMarkers(seurat.all.data, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5)

write.table(seurat.all.data.markers,file="seurat.all.markers.wilcox.res2.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.markers.wilcox.res2.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



head(seurat.all.data.markers)

seurat.all.data.markers %>% group_by(cluster) %>% top_n(10, p_val_adj) -> top10


DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
ggsave("Exp_Heatmap_res2_redblu.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_res2_redblu.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE) + theme(text = element_text(size = 8)) + scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))
ggsave("Exp_Heatmap_res2_BluYR.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_res2_BluYR.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



head(seurat.all.data@meta.data)

Idents(object=seurat.all.data) <- "seurat_clusters"
DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5)  


# assign new cell IDs and stash ID (note: updated with consensus ID)
Idents(object=seurat.all.data) <- "seurat_clusters"
seurat.all.data <- RenameIdents(seurat.all.data, 
`0` = "T cells 1", `1` = "ribosomal-enriched T cells", `2` = "T cells 2", 
`3` = "monocyte/macrophage 1", `4` = "T cells 3", `5` = "activated T cells", 
                                `6` = "monocyte/macrophage 2", 
`7` = "T cells 4", `8` = "heat shock response T cells", `9` = "NK cells 1", 
`10` = "CD8 T cells", `11` = "NK cells 2", `12` = "T cells 5", 
`13` = "T cells 6", `14` = "NK cells 3", `15` = "monocyte/macrophage 3", 
`16` = "XIST+ T cells", `17` = "macrophage 1", `18` = "B cells 1", 
`19` = "T cells 7", `20` = "CD1C+ DCs", `21` = "macrophage 2", 
`22` = "macrophage 3", `23` = "macrophage 4", `24` = "neurons", 
`25` = "CD16 non-classical monocytes", `26` = "proliferating T cells", 
                                `27` = "monocyte/macrophage 4",
`28` = "monocyte/macrophage 5", `29` = "granulocytes", `30` = "monocyte/macrophage 6",
`31` = "B cells 2", `32` = "monocyte-RBCs", `33` = "B cells 3")

seurat.all.data <- StashIdent(object = seurat.all.data, save.name = 'idents_res2_consensus')

head(seurat.all.data@meta.data)

# frequency of cell types across arrays 

#T3
cols.celltypes <- c("#59bd78", # T cells 1
"#59bd78", # ribosomal-enriched T cells
"#59bd78", # T cells 2
"slateblue1", # monocyte/macrophage 1
"#59bd78", # T cells 3
"#59bd78", # activated T cells
"slateblue1", # monocyte/macrophage 2
"#59bd78", # T cells 4
"#59bd78", # heat shock response T cells
"maroon4", # NK cells 1
"yellow2", # CD8 T cells
"maroon4", # NK cells 2
"#59bd78", # T cells 5
"#59bd78", # T cells 6
"maroon4", # NK cells 3
"slateblue1", # monocyte/macrophage 3
"#59bd78", # XIST+ T cells
"#af4add", # macrophage 1
"navy", # B cells 1
"#59bd78", # T cells 7
"turquoise3", # CD1C+ DC
"#af4add", # macrophage 2
"#af4add", # macrophage 3
"#af4add", # macrophage 4
"tomato", # neurons
"slateblue1", # CD16 non-classical monocytes (collapse to monocyte/macrophage)
"#59bd78", # proliferating T cells
"slateblue1", # monocyte/macrophage 4
"slateblue1", # monocyte/macrophage 5
"navajowhite3", # granulocytes
"slateblue1", # monocyte/macrophage 6
"navy", # B cells 2
"slateblue1", # monocyte-RBCs (collapseto monocyte/macrophage)
"navy") # B cells 3


#cols.celltypes <- c("#59bd78", # T cells 1
#"#263e1a", # ribosomal-enriched T cells
#"#59bd78", # T cells 2
#"slateblue1", # monocyte/macrophage 1
#"#59bd78", # T cells 3
#"#6cdd48", # activated T cells
#"slateblue1", # monocyte/macrophage 2
#"#59bd78", # T cells 4
#"#576f47", # heat shock response T cells
#"wheat4", # NK cells 1
#"yellow2", # CD8 T cells
#"wheat4", # NK cells 2
#"#59bd78", # T cells 5
#"#59bd78", # T cells 6
#"wheat4", # NK cells 3
#"slateblue1", # monocyte/macrophage 3
#"#8fab7c", # XIST+ T cells
#"#af4add", # macrophage 1
#"navy", # B cells 1
#"#59bd78", # T cells 7
#"turquoise3", # CD1C+ DC
#"#af4add", # macrophage 2
#"#af4add", # macrophage 3
#"#af4add", # macrophage 4
#"tomato", # neurons
#"#c1b1e1", # CD16 non-classical monocytes
#"#6bce75", # proliferating T cells
#"slateblue1", # monocyte/macrophage 4
#"slateblue1", # monocyte/macrophage 5
#"navajowhite3", # granulocytes
#"slateblue1", # monocyte/macrophage 6
#"navy", # B cells 2
#"#5c3acd", # monocyte-RBCs
#"navy") # B cells 3


# tSNE by new cell types and colors 
#cols.celltypes = c("#5fadb9",
#"#4fbb33",
#"#417171",
#"#880075",
#"#8eb430",
#"#4cb09e",
#"#c640db",
#"#899f3f",
#"#3a7d52",
#"violetred",
#"tomato",
#"violetred2",
#"#3f832d",
#"#8aa78f",
#"violetred4",
#"#765488",
#"#5ab57d",
#"#6c43db",
#"#b29aff",
#"#56682b",
#"turquoise3",
#"#b691c9",
#"#7f6ceb",
#"#7475b6",
#"orangered",
#"#d068dc",
#"yellow3",
#"#8e37b0",
#"#607fe2",
#"navajowhite3",
#"#894d9d",
#"#3c7100",
#"#939ee8",
#"navy")

# all together
Idents(object=seurat.all.data) <- "idents_res2_consensus"

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.celltypes)  
ggsave("tsne_clusterID_annotated.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.celltypes) + NoLegend()
ggsave("tsne_clusterID_annotated_nolegend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_nolegend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



## stacked barplot of clusters by donor, compartment, platform, and hours 
# donor (replicate)
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Replicate)) + geom_bar(position = "fill")
ggsave("StackedBarClusters__res2_replicate.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters__res2_replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# compartment
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Compartment)) + geom_bar(position = "fill")
ggsave("StackedBarClusters__res2_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters__res2_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# platform
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Platform)) + geom_bar(position = "fill")
ggsave("StackedBarClusters__res2_platform.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters__res2_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# hours
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Hours)) + geom_bar(position = "fill")
ggsave("StackedBarClusters__res2_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters__res2_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# array 
ggplot(seurat.all.data@meta.data, aes(x=seurat_clusters, fill=Hours)) + geom_bar(position = "fill")
ggsave("StackedBarClusters__res2_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters__res2_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# proper colors over array 
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + scale_fill_manual(values=cols.celltypes)
ggsave("StackedBarClusters_CellID_Array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_CellID_Array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)




# proper colors over array 
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + scale_fill_manual(values=cols.celltypes) + NoLegend()
ggsave("StackedBarClusters_CellID_Array_nolegend.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_CellID_Array_nolegend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# Separate by compartment (note: need to create correct color palettes for each)
Idents(object=seurat.all.data) <- "Compartment"
suerat.all.data.subset.BLD <- subset(seurat.all.data, idents = c("BLD"))
suerat.all.data.subset.HEF <- subset(seurat.all.data, idents = c("HEF"))


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# delete: macrophage 3, neurons
cols.BLD <- c("#59bd78", # T cells 1
"#59bd78", # ribosomal-enriched T cells
"#59bd78", # T cells 2
"slateblue1", # monocyte/macrophage 1
"#59bd78", # T cells 3
"#59bd78", # activated T cells
"slateblue1", # monocyte/macrophage 2
"#59bd78", # T cells 4
"#59bd78", # heat shock response T cells
"maroon4", # NK cells 1
"yellow2", # CD8 T cells
"maroon4", # NK cells 2
"#59bd78", # T cells 5
"#59bd78", # T cells 6
"maroon4", # NK cells 3
"slateblue1", # monocyte/macrophage 3
"#59bd78", # XIST+ T cells
"#af4add", # macrophage 1
"navy", # B cells 1
"#59bd78", # T cells 7
"turquoise3", # CD1C+ DC
"#af4add", # macrophage 2
"#af4add", # macrophage 4
"slateblue1", # CD16 non-classical monocytes (collapse to monocyte/macrophage)
"#59bd78", # proliferating T cells
"slateblue1", # monocyte/macrophage 4
"slateblue1", # monocyte/macrophage 5
"navajowhite3", # granulocytes
"slateblue1", # monocyte/macrophage 6
"navy", # B cells 2
"slateblue1", # monocyte-RBCs (collapseto monocyte/macrophage)
"navy") # B cells 3


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.BLD)
ggsave("StackedBar_hours_celltypefill_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_hours_celltypefill_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=idents_res2_consensus)) + geom_bar(position = "fill") 


# delete monocyte/macrophage 2, monocyte-RBCs, monocyte/macrophage 5, monocyte/macrophage 4

cols.HF <- c("#59bd78", # T cells 1
"#59bd78", # ribosomal-enriched T cells
"#59bd78", # T cells 2
"slateblue1", # monocyte/macrophage 1
"#59bd78", # T cells 3
"#59bd78", # activated T cells
"#59bd78", # T cells 4
"#59bd78", # heat shock response T cells
"maroon4", # NK cells 1
"yellow2", # CD8 T cells
"maroon4", # NK cells 2
"#59bd78", # T cells 5
"#59bd78", # T cells 6
"maroon4", # NK cells 3
"slateblue1", # monocyte/macrophage 3
"#59bd78", # XIST+ T cells
"#af4add", # macrophage 1
"navy", # B cells 1
"#59bd78", # T cells 7
"turquoise3", # CD1C+ DC
"#af4add", # macrophage 2
"#af4add", # macrophage 3
"#af4add", # macrophage 4
"tomato", # neurons
"slateblue1", # CD16 non-classical monocytes (collapse to monocyte/macrophage)
"#59bd78", # proliferating T cells
"navajowhite3", # granulocytes
"slateblue1", # monocyte/macrophage 6
"navy", # B cells 2
"navy") # B cells 3


ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=idents_res2_consensus)) + geom_bar(position = "fill") + scale_fill_manual(values=cols.HF)
ggsave("StackedBar_hours_celltypefill_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_hours_celltypefill_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.BLD)
ggsave("StackedBar_arrayname_celltypefill_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_celltypefill_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=ArrayName, fill=idents_res2_consensus)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=cols.HF)
ggsave("StackedBar_arrayname_celltypefill_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_celltypefill_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


seurat.all.factors <- seurat.all.data@meta.data$ArrayName
plate.table = table(seurat.all.factors, as.vector(seurat.all.data@meta.data$idents_res2_collapse_final))
head(plate.table)
write.table(plate.table,file="cellfreqs.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp cellfreqs.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_alldata_filtered.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


help(VlnPlot)

# violin plots of top marker genes 
Idents(object=seurat.all.data) <- "idents_res2_consensus"
VlnPlot(seurat.all.data, features = c("PDCD4","RPS27","IL32","IL1B",
                                     "RPSAP58","ICOS","S100A9","MTRNR2L1"), pt.size = 0, y.max = 8, cols = cols.celltypes)
ggsave("Marker_Violin_ClusterID.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Marker_Violin_ClusterID.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c("HSPA1A","PRF1","CCL5","PTMA","CCL2","TIPIN","FGFBP2",
                                      "HLA-DRA"), pt.size = 0, y.max = 8, cols = cols.celltypes)
ggsave("Marker_Violin_ClusterID2.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Marker_Violin_ClusterID2.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


VlnPlot(seurat.all.data, features = c("MALAT1","THBS1","IGHM","ELMSAN1","FCER1A",
                                     "CTSL","CXCL5","IL8"), pt.size = 0, y.max = 8, cols = cols.celltypes)
ggsave("Marker_Violin_ClusterID3.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Marker_Violin_ClusterID3.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("KIF5A","LST1","STMN1","IFI44L","FCN1","CLC","FGFBP2","IGJ","GPM6A","IGHG2"), pt.size = 0, y.max = 8, cols = cols.celltypes)
ggsave("Marker_Violin_ClusterID4.pdf", useDingbats = F, height = 12, width = 16, units = "in")
system(paste0("gsutil cp Marker_Violin_ClusterID4.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# QC violins by cell type 
Idents(object=seurat.all.data) <- "idents_res2_consensus"

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1, cols = cols.celltypes)
ggsave("nFeatureRNA_violin_annotated.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nFeatureRNA_violin_annotated.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1,cols = cols.celltypes)
ggsave("nCountRNA_violin_annotated.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nCountRNA_violin_annotated.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1,cols = cols.celltypes)
ggsave("ptMT_violin_annotated.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp ptMT_violin_annotated.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0, cols = cols.celltypes)
ggsave("nFeatureRNA_violin_annotated_nopoints.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nFeatureRNA_violin_annotated_nopoints.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0,cols = cols.celltypes)
ggsave("nCountRNA_violin_annotated_nopoints.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nCountRNA_violin_annotated_nopoints.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0,cols = cols.celltypes)
ggsave("ptMT_violin_annotated_nopoints.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp ptMT_violin_annotated_nopoints.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)




head(seurat.all.data@meta.data)

# assign new cell IDs and stash ID (note: updated with collapsed ID - also collapsed macs and monos)
#Idents(object=seurat.all.data) <- "RNA_snn_res.2"
#seurat.all.data <- RenameIdents(seurat.all.data, 
#`0` = "T cells", `1` = "T cells", `2` = "T cells", 
#`3` = "monocyte/macrophage", `4` = "T cells", `5` = "T cells", 
#                                `6` = "monocyte/macrophage", 
#`7` = "T cells", `8` = "T cells", `9` = "NK cells", 
#`10` = "CD8 T cells", `11` = "NK cells", `12` = "T cells", 
#`13` = "T cells", `14` = "NK cells", `15` = "monocyte/macrophage", 
#`16` = "T cells", `17` = "macrophage", `18` = "B cells", 
#`19` = "T cells", `20` = "CD1C+ DCs", `21` = "macrophage", 
#`22` = "macrophage", `23` = "macrophage", `24` = "neurons", 
#`25` = "monocyte/macrophage", `26` = "T cells", 
#                                `27` = "monocyte/macrophage",
#`28` = "monocyte/macrophage", `29` = "granulocytes", `30` = "monocyte/macrophage",
#`31` = "B cells", `32` = "monocyte/macrophage", `33` = "B cells")

# assign new cell IDs and stash ID (note: updated with collapsed ID - also collapsed macs and monos)
Idents(object=seurat.all.data) <- "RNA_snn_res.2"
seurat.all.data <- RenameIdents(seurat.all.data, 
`0` = "T cells", `1` = "T cells", `2` = "T cells", 
`3` = "monocyte", `4` = "T cells", `5` = "T cells", 
`6` = "monocyte", `7` = "T cells", `8` = "T cells", 
`9` = "NK cells", `10` = "CD8 T cells", `11` = "NK cells", 
`12` = "T cells", `13` = "T cells", `14` = "NK cells", `15` = "monocyte", 
`16` = "T cells", `17` = "monocyte", `18` = "B cells", 
`19` = "T cells", `20` = "CD1C+ DCs", `21` = "monocyte", 
`22` = "monocyte", `23` = "monocyte", `24` = "neurons", 
`25` = "monocyte", `26` = "T cells", 
`27` = "monocyte", `28` = "monocyte", `29` = "granulocytes", 
`30` = "monocyte",
`31` = "B cells", `32` = "monocyte", `33` = "B cells")

#seurat.all.data <- StashIdent(object = seurat.all.data, save.name = 'idents_res2_colapse')
seurat.all.data <- StashIdent(object = seurat.all.data, save.name = 'idents_res2_collapse_final')

unique(seurat.all.data@meta.data$idents_res2_collapse_final)

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1)


cols.celltypes.collapse <- c("#59bd78", # T cells
                             "slateblue1", # monocyte
                             "maroon4", # NK cells
                             "yellow2", # CD8 T cells
                            "navy", # B cells
                            "turquoise3", # CD1C+ DC
                            "tomato", # neurons
                            "navajowhite3") # granulocytes  
#"#af4add", # macrophage 2

Idents(object=seurat.all.data) <- "idents_res2_collapse_final"


unique(seurat.all.data@meta.data$idents_res2_collapse_final)

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0.1, cols = cols.celltypes.collapse)
ggsave("nFeatureRNA_violin_annotated_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nFeatureRNA_violin_annotated_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0.1,cols = cols.celltypes.collapse)
ggsave("nCountRNA_violin_annotated_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nCountRNA_violin_annotated_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0.1,cols = cols.celltypes.collapse)
ggsave("ptMT_violin_annotated_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp ptMT_violin_annotated_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nFeature_RNA"), pt.size = 0, cols = cols.celltypes.collapse)
ggsave("nFeatureRNA_violin_annotated_nopoints_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nFeatureRNA_violin_annotated_nopoints_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("nCount_RNA"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("nCountRNA_violin_annotated_nopoints_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp nCountRNA_violin_annotated_nopoints_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("ptMT_violin_annotated_nopoints_colapse.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp ptMT_violin_annotated_nopoints_colapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# tSNE colored by compartment, separated by time
DimPlot(seurat.all.data, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Hours")
ggsave("UMAP_hours_splitbyhrs_fixed.pdf", useDingbats = F, height = 8, width = 16, units = "in")
system(paste0("gsutil cp UMAP_hours_splitbyhrs_fixed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# tSNE colored by compartment, separated by time
DimPlot(suerat.all.data.subset.BLD, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Hours")
ggsave("tsne_splitbyhrs_BLD.pdf", useDingbats = F, height = 8, width = 16, units = "in")
system(paste0("gsutil cp tsne_splitbyhrs_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# tSNE colored by compartment, separated by time
DimPlot(suerat.all.data.subset.HEF, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Hours")
ggsave("tsne_splitbyhrs_HEF.pdf", useDingbats = F, height = 8, width = 16, units = "in")
system(paste0("gsutil cp tsne_splitbyhrs_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


head(suerat.all.data.subset.BLD@meta.data)

# tSNE colored by compartment, separated by time
DimPlot(suerat.all.data.subset.BLD, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Hours")
ggsave("tsne_splitbyarray_BLD.pdf", useDingbats = F, height = 8, width = 32, units = "in")
system(paste0("gsutil cp tsne_splitbyarray_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# Separate by compartment (note: need to create correct color palettes for each)
Idents(object=suerat.all.data.subset.BLD) <- "Hours"
suerat.all.data.subset.BLD.PT <- subset(suerat.all.data.subset.BLD, idents = c("050","066","073","089","112","137"))


# tSNE colored by compartment, separated by time
DimPlot(suerat.all.data.subset.BLD.PT, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = cols.hours, split.by = "Hours")
ggsave("tsne_splitbyarray_BLD_NF.pdf", useDingbats = F, height = 8, width = 16, units = "in")
system(paste0("gsutil cp tsne_splitbyarray_BLD_NF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


Idents(object=suerat.all.data.subset.BLD) <- "Hours"
suerat.all.data.subset.BLD.follow <- subset(suerat.all.data.subset.BLD, idents = c("500"))


# tSNE colored by compartment, separated by time
DimPlot(suerat.all.data.subset.BLD.follow, reduction = "tsne", group.by = "Hours", pt.size = 0.8, cols = c("gold"), split.by = "ArrayName")
ggsave("tsne_splitbyarray_BLD_follow.pdf", useDingbats = F, height = 8, width = 16, units = "in")
system(paste0("gsutil cp tsne_splitbyarray_BLD_follow.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_alldata_filtered.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# save blood to analyze seperately  
saveRDS(suerat.all.data.subset.BLD,file="seuratobj_alldata_filtered_BLD.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_BLD.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# save hematoma to analyze seperately
saveRDS(suerat.all.data.subset.HEF,file="seuratobj_alldata_filtered_HEF.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_HEF.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


unique(seurat.all.data@meta.data$idents_res2_consensus)

# subset out T cells 
Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.subset.Tcells <- subset(seurat.all.data, idents = c("activated T cells",
"proliferating T cells",
"CD8 T cells",
"NK cells 1",
"heat shock response T cells",
"T cells 4",
"NK cells 3",
"NK cells 2",
"T cells 2",
"ribosomal-enriched T cells",
"XIST+ T cells",
"T cells 1",
"T cells 5",
"T cells 6",
"T cells 7"))
saveRDS(suerat.all.data.subset.Tcells,file="seuratobj_alldata_filtered_Tcells.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_Tcells.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# subset out myeloid cells
Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.subset.myeloid <- subset(seurat.all.data, idents = c("monocyte/macrophage 1",
"CD16 non-classical monocytes",
"CD1C+ DCs",
"monocyte/macrophage 6",
"granulocytes",
"monocyte/macrophage 2",
"macrophage 4",
"macrophage 1",
"monocyte/macrophage 3",
"monocyte/macrophage 4",
"monocyte-RBCs",
"monocyte/macrophage 5",
"macrophage 3",
"macrophage 2"))
saveRDS(suerat.all.data.subset.myeloid,file="seuratobj_alldata_filtered_myeloid.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_myeloid.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# NK cells
Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.subset.NK <- subset(seurat.all.data, idents = c("NK cells 1","NK cells 2","NK cells 3"))
saveRDS(suerat.all.data.subset.NK,file="seuratobj_alldata_filtered_NK.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_NK.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


# neurons 
Idents(object=seurat.all.data) <- "idents_res2_consensus"
suerat.all.data.subset.neurons <- subset(seurat.all.data, idents = c("neurons"))
saveRDS(suerat.all.data.subset.neurons,file="seuratobj_alldata_filtered_neurons.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered_neurons.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


Idents(object=seurat.all.data) <- "idents_res2_collapse_final"
VlnPlot(seurat.all.data, features = c("percent.mt"), pt.size = 0,cols = cols.celltypes.collapse)


Idents(object=seurat.all.data) <- "idents_res2_collapse_final"
DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.celltypes.collapse) + NoLegend()
ggsave("tsne_clusterID_annotated_collapse_nolegend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_collapse_nolegend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.celltypes.collapse)  
ggsave("tsne_clusterID_annotated_collapse.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_collapse.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.celltypes.collapse) + NoLegend()
ggsave("tsne_clusterID_annotated_collapse_nolegend.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_collapse_nolegend.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


## stacked barplot of clusters by donor, compartment, platform, and hours 
# donor (replicate)
ggplot(seurat.all.data@meta.data, aes(x=idents_res2_collapse_final, fill=Replicate)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90))
ggsave("StackedBarClusters_collapse_replicate.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_replicate.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# compartment
ggplot(seurat.all.data@meta.data, aes(x=idents_res2_collapse_final, fill=Compartment)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90)) + scale_fill_manual(values=cols.compartment)
ggsave("StackedBarClusters_collapse_compartment.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_compartment.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# platform
ggplot(seurat.all.data@meta.data, aes(x=idents_res2_collapse_final, fill=Platform)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90)) + scale_fill_manual(values=cols.platform)
ggsave("StackedBarClusters_collapse_platform.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_platform.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# hours
ggplot(seurat.all.data@meta.data, aes(x=idents_res2_collapse_final, fill=Hours)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90)) + scale_fill_manual(values=cols.hours) 
ggsave("StackedBarClusters_collapse_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# array 
ggplot(seurat.all.data@meta.data, aes(x=idents_res2_collapse_final, fill=ArrayName)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90))
ggsave("StackedBarClusters_collapse_Hours.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_Hours.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# proper colors over array 
ggplot(seurat.all.data@meta.data, aes(x=ArrayName, fill=idents_res2_collapse_final)) + geom_bar(position = "fill") + scale_fill_manual(values=cols.celltypes.collapse) + theme(axis.text.x = element_text(angle = 90))
ggsave("StackedBarClusters_collapse_array.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBarClusters_collapse_array.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=idents_res2_collapse_final)) + 
geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

cols.celltypes.collapse.BLD <- c("#59bd78", # T cells
                             "slateblue1", # monocyte
                             "maroon4", # NK cells 1
                             "yellow2", # CD8 T cells
                            "navy", # B cells 2
                            "turquoise3", # CD1C+ DC
                            "navajowhite3") # granulocytes 

ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=Hours, fill=idents_res2_collapse_final)) + 
geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_fill_manual(values = cols.celltypes.collapse.BLD)
ggsave("StackedBar_hours_collapse_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_hours_collapse_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=Hours, fill=idents_res2_collapse_final)) + 
geom_bar(position = "fill") + scale_fill_manual(values = cols.celltypes.collapse)
ggsave("StackedBar_hours_collapse_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_hours_collapse_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.BLD@meta.data, aes(x=ArrayName, fill=idents_res2_collapse_final)) + geom_bar(position = "fill") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values=cols.celltypes.collapse.BLD)
ggsave("StackedBar_arrayname_collapse_BLD.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_collapse_BLD.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.HEF@meta.data, aes(x=ArrayName, fill=idents_res2_collapse_final)) + 
geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values=cols.celltypes.collapse)
ggsave("StackedBar_arrayname_collapse_HEF.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_collapse_HEF.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


saveRDS(seurat.all.data,file="seuratobj_alldata_filtered.rds")
system(paste0("gsutil cp seuratobj_alldata_filtered.rds gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


head(seurat.all.data@meta.data)

Idents(object=seurat.all.data) <- "idents_res2_collapse_final"

seurat.all.data.markers <- FindAllMarkers(seurat.all.data, test.use = "wilcox", 
                                          only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5)


write.table(seurat.all.data.markers,file="seurat.all.markers.wilcox.collapsedclusters.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp seurat.all.markers.wilcox.collapsedclusters.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



#seurat.all.data.markers %>% group_by(cluster) %>% top_n(10, p_val_adj) -> top10

top10 <- seurat.all.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


head(top10)

DoHeatmap(seurat.all.data,features=top10$gene, raster = FALSE)
ggsave("Exp_Heatmap_CollapsedCellIDs.pdf", useDingbats = F)
system(paste0("gsutil cp Exp_Heatmap_CollapsedCellIDs.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



# T cells
VlnPlot(seurat.all.data, features = c("IL7R"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_Tcells_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_Tcells_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# Monocyte
VlnPlot(seurat.all.data, features = c("S100A9"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_MM_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_MM_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# NK cells
VlnPlot(seurat.all.data, features = c("GNLY"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_NKcells_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_NKcells_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# CD8 T cells
VlnPlot(seurat.all.data, features = c("GZMK"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_CD8Tcells_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_CD8Tcells_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#B cells
VlnPlot(seurat.all.data, features = c("IGKC"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_Bcells_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_Bcells_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# DCs
VlnPlot(seurat.all.data, features = c("HLA-DRA"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_CD1CDCs_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_CD1CDCs_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

# neurons
VlnPlot(seurat.all.data, features = c("KIF5A"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_Neurons_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_Neurons_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

#grans
VlnPlot(seurat.all.data, features = c("CLC"), pt.size = 0,cols = cols.celltypes.collapse)
ggsave("Violin_grans_Collapsed.pdf", useDingbats = F, height = 7, width = 14, units = "in")
system(paste0("gsutil cp Violin_grans_Collapsed.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


#Idents(object=seurat.all.data) <- "ArrayName"
seurat.all.factors.collapse <- seurat.all.data@meta.data$ArrayName
plate.table.collapse = table(seurat.all.factors.collapse, as.vector(seurat.all.data@meta.data$idents_res2_colapse))
head(plate.table.collapse)
write.table(plate.table.collapse,file="cellfreqs_collapsed.txt",sep="\t",quote=FALSE,col.names=NA)
system(paste0("gsutil cp cellfreqs_collapsed.txt gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


unique(seurat.all.data@meta.data$ArrayName)

# subset out just control blood 

# Separate by compartment (note: need to create correct color palettes for each)
Idents(object=seurat.all.data) <- "ArrayName"
suerat.all.data.subset.control <- subset(seurat.all.data, idents = c("DayF-500-BLD-00-Nova","DayF-500-BLD-01-Nova",
                                                                 "DayF-500-BLD-02-Nova","DayF-500-BLD-03-Nova"))


cols.celltypes.collapse.controls <- c("#59bd78", # T cells
                             "slateblue1", # monocyte/macrophage 1
                             "maroon4", # NK cells 1
                             "yellow2", # CD8 T cells
                            "navy", # B cells 2
                            "turquoise3", # CD1C+ DC
                            "navajowhite3") # granulocytes                         

ggplot(suerat.all.data.subset.control@meta.data, aes(x=ArrayName, fill=idents_res2_collapse_final)) + geom_bar(position = "fill") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values=cols.celltypes.collapse.controls)
ggsave("StackedBar_arrayname_collapse_controls.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_collapse_controls.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


ggplot(suerat.all.data.subset.control@meta.data, aes(x=Hours, fill=idents_res2_collapse_final)) + geom_bar(position = "fill") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values=cols.celltypes.collapse.controls)
ggsave("StackedBar_hours_collapse_controls.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_hours_collapse_controls.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


Idents(object=seurat.all.data) <- "ArrayName"
suerat.all.data.subset.BLD.CT <- subset(seurat.all.data, idents = c("DayF-500-BLD-PT-Nova"))

ggplot(suerat.all.data.subset.BLD.CT@meta.data, aes(x=ArrayName, fill=idents_res2_collapse_final)) + geom_bar(position = "fill") + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
scale_fill_manual(values=cols.celltypes.collapse.controls)
ggsave("StackedBar_arrayname_collapse_PT.pdf", useDingbats = F)
system(paste0("gsutil cp StackedBar_arrayname_collapse_PT.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)


Idents(object=seurat.all.data) <- "ArrayName"

unique(seurat.all.data@meta.data$ArrayName)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5)  


#cols.hours = c(05 = "lightsteelblue1", 66 = "darkslategray1",
# 73 = "cadetblue2", 89 = "dodgerblue2",112 = "blue1", 137 = "darkblue",controls = "gold", PT = "orange")

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
ggsave("tsne_clusterID_annotated_collapse_PThighlighted.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_collapse_PThighlighted.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)

DimPlot(seurat.all.data, reduction = "tsne", pt.size = 0.5, cols = cols.array.pt.ctl) + NoLegend()
ggsave("tsne_clusterID_annotated_collapse_nolegend_PThighlighted.pdf", useDingbats = F, height = 8, width = 8, units = "in")
system(paste0("gsutil cp tsne_clusterID_annotated_collapse_nolegend_PThighlighted.pdf gs://fc-a6329ed4-8283-48d4-90c3-d7175f18a168/Analysis/20200401_MultipletsRemoved_Figure1/"), intern=TRUE)



