#FixedRNA_Yan_Feb2024 - Single cell transcriptomic analysis of Neutrophils on COVID19 patients
#March 2024 -Kalpani de Silva

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(gprofiler2)
library(ggrepel)
library(dplyr) 
library(scales)
library(SeuratObject)
library(DropletUtils)
library(SingleCellExperiment)

#setwd("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/GSE205013")
setwd("/bio/home/kkdesi01/KBRIN0XXX-Yan-FixedRNA")

#~~~~~~~~~~~~~~~~~~~~~~
#Merging with EmptyDrops
#~~~~~~~~~~~~~~~~~~~~~~


#"Post_C19-1_CD16_high" --> "Post_C19_1_CD16_high","Post_C19-1_NDN" --> "Post_C19_1_NDN"

filelist<-c("C19_P1_NDN","C19_P2_NDN","C19_P3_NDN","C19_P4_NDN","C19_P5_NDN","C19_P6_NDN","C19_P7_NDN","Healthy_d1_NDN","Healthy_d2_NDN","Healthy_d3_NDN","Healthy_d4_NDN","Healthy_d5_NDN","Non_C19_P1_NDN","Non_C19_P2_NDN")

grouplist<-c("C19_NDN","C19_NDN","C19_NDN","C19_NDN","C19_NDN","C19_NDN","C19_NDN","Healthy_NDN","Healthy_NDN","Healthy_NDN","Healthy_NDN","Healthy_NDN","Non_C19_NDN","Non_C19_NDN")

celltype <-c("NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN","NDN")

#Create seurat object for each sample
for (i in filelist) {
  sce <- Read10X(data.dir = paste0("./CellRanger/", i, "/count/sample_raw_feature_bc_matrix/"))
  set.seed(100)
  limit <- 100
  out <- emptyDrops(sce, lower=limit)
  sce <- sce[,which(out$FDR <= 0.001)]

  covid.seurat.obj <- CreateSeuratObject(counts = sce, project = "covid")
  covid.seurat.obj@meta.data$sample <- i
  covid.seurat.obj@meta.data$group <- grouplist[match(i, filelist)]
  covid.seurat.obj@meta.data$celltype <- celltype[match(i, filelist)]
  assign(i, covid.seurat.obj)
}

print("Seurat objects generated for each sample")


#Merge all Seurat objects in the list
covid.seurat.obj <-  merge(C19_P1_NDN, y = c(C19_P2_NDN,C19_P3_NDN,C19_P4_NDN,C19_P5_NDN,C19_P6_NDN,C19_P7_NDN,Healthy_d1_NDN,Healthy_d2_NDN,Healthy_d3_NDN,Healthy_d4_NDN,Healthy_d5_NDN,Non_C19_P1_NDN,Non_C19_P2_NDN),
                           add.cell.ids = c("C19_P1_NDN","C19_P2_NDN","C19_P3_NDN","C19_P4_NDN","C19_P5_NDN","C19_P6_NDN","C19_P7_NDN","Healthy_d1_NDN","Healthy_d2_NDN","Healthy_d3_NDN","Healthy_d4_NDN","Healthy_d5_NDN","Non_C19_P1_NDN","Non_C19_P2_NDN"), project = "COVID19")


#sanity check to see whether we have all the patients data merged
unique(covid.seurat.obj@meta.data$sample)

unique(covid.seurat.obj@meta.data$group)

unique(covid.seurat.obj@meta.data$celltype)

#how many cells left after Emptydrop filtering
sample_counts <- table(covid.seurat.obj$sample)
print(sample_counts)

# Specify the correct file path and name for saving the .RDS file
file_path <- "./NDN_Analysis/my_rds/covid.seurat.obj.EmptyDrops.rds"

# Save the Seurat object as .RDS file
saveRDS(covid.seurat.obj, file = file_path)

covid.seurat.obj

print("Seurat objects are merged into covid.seurat.obj")


#~~~~~~~~~~~~~~~~~~~~~~
#Filtering
#~~~~~~~~~~~~~~~~~~~~~~
covid.seurat.obj <- readRDS("./NDN_Analysis/my_rds/covid.seurat.obj.EmptyDrops.rds")

# calculate mitochondrial gene percentage
covid.seurat.obj$MTPercent <- PercentageFeatureSet(covid.seurat.obj, pattern = '^MT-' )

##some columns have NaN for MTPercent. Cant create vln plot with NaN.
## So I am transforming NaN into 0

covid.seurat.obj@meta.data$MTPercent[is.nan(covid.seurat.obj@meta.data$MTPercent)] <- 0

# Example: Impute missing values with median
#pdac.seurat.obj$MTPercent[is.na(pdac.seurat.obj$MTPercent)] <- median(pdac.seurat.obj$MTPercent, na.rm = TRUE)
# Example: Remove cells with missing MT percentage
#pdac.seurat.obj <- subset(pdac.seurat.obj, subset = !is.na(MTPercent))

# calculate ribosomal gene percentage
covid.seurat.obj$RiboPercent <- PercentageFeatureSet(covid.seurat.obj, pattern = '^RP[SL]')
covid.seurat.obj@meta.data$RiboPercent[is.nan(covid.seurat.obj@meta.data$RiboPercent)] <- 0
# Example: Impute missing values with median
#pdac.seurat.obj$RiboPercent[is.na(pdac.seurat.obj$RiboPercent)] <- median(pdac.seurat.obj$RiboPercent, na.rm = TRUE)

#calculating CD66b gene percentages inorder to keep only the neutrophils
#genes <- c("CEACAM8", "ITGAM", "FCGR3A", "FUT4", "SELL", "CD177", "CD63", "ELANE")\
#genes <- c("MS4A1") #CD20
#genes <- c("CD3")
genes <- c("MS4A1", "CD3D") #"MS4A1"- B cells, "CD3D"- T cells

covid.seurat.obj$MarkerPercent <- PercentageFeatureSet(covid.seurat.obj, features = genes )

#THIS IS TOO BIG CANNOT CREATE
#output_file1 <- "./New_Analysis/my_rds/1Filtering/Violinplot_beforeFiltering_AfterEmptyDrops_CD66b.pdf"
#output_format <- "pdf"
#ggsave(filename = output_file1, plot = violinplot1, device = output_format)

#featurescatter <- FeatureScatter(covid.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE) + geom_smooth(method = 'lm')
#output_file2 <- "./New_Analysis/my_rds/1Filtering/featurescatter_beforeFiltering_AfterEmptyDrops.pdf"
#output_format <- "pdf"
#ggsave(filename = output_file2, plot = featurescatter, device = output_format)

# filtering

covid.seurat.obj <- subset(covid.seurat.obj, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & MTPercent < 10 & MarkerPercent <= 0 )

violinplot2 <- VlnPlot(covid.seurat.obj, features = c("nCount_RNA", "nFeature_RNA", "MTPercent", "MarkerPercent"), ncol = 4)
output_file3 <- "./NDN_Analysis/my_rds/1Filtering/violinplot_AfterFiltering_AfterEmptyDrops_Bfilter.pdf"
output_format <- "pdf"
ggsave(filename = output_file3, plot = violinplot2, device = output_format)

#Normalization
covid.seurat.obj <- NormalizeData(covid.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify highly variable features
covid.seurat.obj <- FindVariableFeatures(covid.seurat.obj, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(covid.seurat.obj), 20)
plot1 <- VariableFeaturePlot(covid.seurat.obj)
plot1 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)


output_file4 <- "./NDN_Analysis/my_rds/1Filtering/VariableFeaturesAfterFiltering_AfterEmptyDrops_Bfilter.pdf"
output_format <- "pdf"
ggsave(filename = output_file4, plot = plot1, device = output_format)



covid.seurat.obj


# Specify the correct file path and name for saving the .RDS file to move to cluster for integration step
file_path <- "./NDN_Analysis/my_rds/covid.filtered.Bfilter.seurat.obj.EmptyDrops.rds"

#Save the Seurat object as .RDS file
saveRDS(covid.seurat.obj, file = file_path)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#After saving the covid.filtered.seurat.obj.rds, cluster it and generate a UMAP to see if clusters are bias to patients

#Scaling
covid.filtered.seurat.obj <- ScaleData(object = covid.seurat.obj)

#Perform Linear dimensionality reduction
covid.filtered.seurat.obj  <- RunPCA(object = covid.filtered.seurat.obj )
# visualize PCA results
print(covid.filtered.seurat.obj [["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(covid.filtered.seurat.obj , dims = 1:5, cells = 500, balanced = TRUE)
#selecting statistically significant PCs
Elbowplot <- ElbowPlot(covid.filtered.seurat.obj , ndims = 40)
output_file10 <- "./NDN_Analysis/my_rds/1Filtering/ElbowplotAfterFilteringBeforeIntegration_EmptyDrops_TBfilter.pdf"
output_format <- "pdf"
ggsave(filename = output_file10, plot = Elbowplot, device = output_format)


#Clustering
covid.filtered.seurat.obj <- FindNeighbors(covid.filtered.seurat.obj, dims = 1:30)
covid.filtered.seurat.obj <- FindClusters(covid.filtered.seurat.obj, resolution = c(0.1))


#my_palette <- c('#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#FFA500', '#808080', '#FFC0CB')

# Assign identity of clusters
Idents(object = covid.filtered.seurat.obj) <- "RNA_snn_res.0.1"
covid.filtered.seurat.obj <- RunUMAP(object = covid.filtered.seurat.obj, dims = 1:30)
p3 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'sample')
p31 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'group')
p32 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'celltype')
p33 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
grid3 <- grid.arrange(p3, p31, p32, p33, ncol = 2, nrow = 2)



# Save the plot to the specified file
output_file6 <- "./NDN_Analysis/my_rds/1Filtering/EmptyDrops_UMAPAfterFilteringRes0.1_TBfilter.pdf"
output_format <- "pdf"
pdf_width <- 30
pdf_height <- 25
ggsave(filename = output_file6, plot = grid3, device = output_format, width = pdf_width, height = pdf_height)


#how many cells left after EmptyDrop and Manual filtering
sample_counts <- table(covid.filtered.seurat.obj$sample)
print(sample_counts)

group_counts <- table(covid.filtered.seurat.obj$group)
print(group_counts)

celltype_counts <- table(covid.filtered.seurat.obj$celltype)
print(celltype_counts)

seurat_clusters_counts <- table(covid.filtered.seurat.obj$seurat_clusters)
print(seurat_clusters_counts)



sessionInfo()