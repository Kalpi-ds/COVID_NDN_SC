#FixedRNA_Feb2024 - Single cell transcriptomic analysis of Neutrophils on COVID19 patients
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


setwd("/path/KBRIN0XXX-xxx-FixedRNA")


#~~~~~~~~~~~~~~~~~~~~~~
#FindMarkers
#~~~~~~~~~~~~~~~~~~~~~~

# Load the RDS file
covid.seurat.obj <- readRDS("./NDN_Analysis/my_rds/covid.filtered.seurat.HarmonyIntegration.obj.EmptyDrops.BTfilter.rds")

#use the identified best resolution for clustering
Idents(object = covid.seurat.obj) <- "seurat_clusters"


# cluster re-assignment occurs, which re-assigns clustering in my_levels
#my_levels <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

# Re-level object@ident
#covid.seurat.obj@active.ident <- factor(x = covid.seurat.obj@active.ident, levels = my_levels)


#covid.seurat.obj <- RunUMAP(object = covid.seurat.obj, dims = 1:30)
p1 <- DimPlot(covid.seurat.obj, reduction = 'umap', label = TRUE, label.size = 6)
# Save the plot to the specified file
output_file1 <- "./NDN_Analysis/my_rds/4Markers/UMAP_Res0.3.pdf"
output_format <- "pdf"
ggsave(filename = output_file1, plot = p1, device = output_format)

#change the default assay into RNA
DefaultAssay(covid.seurat.obj)
DefaultAssay(covid.seurat.obj) <- "RNA"
DefaultAssay(covid.seurat.obj)

#scale.data slot for the RNA assay
covid.seurat.obj <- ScaleData(covid.seurat.obj)


#find markers with log2fc 0.25
All_Markers <- FindAllMarkers(covid.seurat.obj, only.pos=FALSE, assay = "RNA")
All_Markers <- All_Markers %>% group_by(cluster)
write.table(All_Markers, file = "./NDN_Analysis/my_rds/4Markers/AllMarkers.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Group by cluster and select the top 10 markers based on avg_log2FC, %>% = pipe"|"
top10 <- All_Markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
# Write the top 10 markers to a tab-delimited text file
write.table(top10, file = "./NDN_Analysis/my_rds/4Markers/Top10AllMarkers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#heatmap of top10 markers
Top10AllMarkers <- DoHeatmap(covid.seurat.obj, features = top10$gene)
# Modify the y-axis gene labels font size and style
Top10AllMarkers <- Top10AllMarkers + theme(axis.text.y = element_text(size = 11, face = "bold"))
pdf_width <- 30
pdf_height <- 45
output_file2 <- "./NDN_Analysis/my_rds/4Markers/Top10AllMarkersHeatmap.pdf"
output_format <- "pdf"
ggsave(filename = output_file2, plot = Top10AllMarkers, device = output_format, width = pdf_width, height = pdf_height)


# Create a function to generate dimplot for each cluster
for (cluster_id in unique(top10$cluster)){
  cluster_markers <- filter(top10, cluster == cluster_id)$gene
  Top10Dimplots <- FeaturePlot(covid.seurat.obj, features = cluster_markers)

  # Set the desired width and height for the PDF
  pdf_width <- 45
  pdf_height <- 45

  output_file3 <- paste("./NDN_Analysis/my_rds/4Markers/Cluster", cluster_id, "Top10Dimplots.pdf", sep = "")
  output_format <- "pdf"
  ggsave(filename = output_file3, plot = Top10Dimplots, device = output_format, width = pdf_width, height = pdf_height)
}

# Create a function to generate violin plots for each cluster's top markers
for (cluster_id in unique(top10$cluster)){
  cluster_markers <- filter(top10, cluster == cluster_id)$gene
  Top10Vlnplots <- VlnPlot(covid.seurat.obj, features = cluster_markers, ncol = 2)

  # Set the desired width and height for the PDF
  pdf_width <- 35  # Modify this value to adjust the width of the PDF in inches
  pdf_height <- 45  # Modify this value to adjust the height of the PDF in inches

  output_file4 <- paste("./NDN_Analysis/my_rds/4Markers/Cluster", cluster_id, "Top10Vlnplots.pdf", sep = "")
  output_format <- "pdf"
  ggsave(filename = output_file4, plot = Top10Vlnplots, device = output_format, width = pdf_width, height = pdf_height)
}



#Get the number of cells in each cluster by Treatment
Table <- table(Idents(object = covid.seurat.obj))
write.table(Table, file = "./NDN_Analysis/my_rds/4Markers/numberofcellsineachclusterbyTreatment.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#find all markers in the clusters
#getting all the markers with out log2fc filtering. so we can see all the genes in the data set. Switch back to RNA assay when doing differential expression
All_Markers_nofilter <- FindAllMarkers(covid.seurat.obj, only.pos=FALSE, logfc.threshold = 0, assay = "RNA")
All_Markers_nofilter <- All_Markers_nofilter %>% group_by(cluster)
write.table(All_Markers_nofilter, file = "./NDN_Analysis/my_rds/4Markers/AllMarkers_nofilters.txt", sep = "\t", quote = FALSE, row.names = FALSE)


covid.seurat.obj

# Save the Seurat object as .RDS file
file_path <- "./NDN_Analysis/my_rds/covid.NDN.0.3.Markers.BTfilter.seurat.obj.rds"

saveRDS(covid.seurat.obj, file = file_path)

