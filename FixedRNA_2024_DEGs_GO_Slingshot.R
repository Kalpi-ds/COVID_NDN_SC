
#FixedRNA_Yan_Feb2024 - Single cell transcriptomic analysis of Neutrophils on COVID19 patients
#June 2024 -Kalpani de Silva


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(gprofiler2)
library(ggrepel)
library(dplyr) 
library(scales)
library(slingshot)


setwd("/bio/home/kkdesi01/KBRIN0XXX-Yan-FixedRNA")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DEGs in NDN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NDN_subset.seurat.obj <- readRDS("./NDN_Analysis/my_rds/covid.NDN.0.3.Markers.BTfilter.seurat.obj.rds")
# 
# NDN_subset.seurat.obj
# 
# DefaultAssay(NDN_subset.seurat.obj) <- "RNA" #The default when running ScaleData is that only the variable features will be scaled. To scale all features all you need to do is set assay to "RNA" and run
# NDN_subset.seurat.obj <- ScaleData(NDN_subset.seurat.obj, features = rownames(NDN_subset.seurat.obj))
# 
# 
# #change to the appropriate cluster number
# Idents(object = NDN_subset.seurat.obj) <- "group"
# 
# 
# #C19_NDN_VS_Non_C19_NDN
# 
# for(i in 0:7){
# 
#   #THE OPTIONS ARE SET for no filtering of results.
#   seurat_object_sub<-subset(NDN_subset.seurat.obj, subset = seurat_clusters == i)
#   markers <- FindMarkers(seurat_object_sub, ident.1 = "C19_NDN",ident.2="Non_C19_NDN", only.pos=FALSE, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, assay = "RNA")
#   fn<-paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster",i,"_NoFiltering.txt",sep="")
#   write.table(markers,fn,sep="\t", col.names = NA)
# 
#   #Filter markers by p_val_adj < 0.05
#   filtered_markers <- markers[markers$p_val_adj < 0.05, ]
# 
#   fn2 <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
#   write.table(filtered_markers, fn2, sep = "\t", col.names = NA)
# }
# 
# 
# for(i in 0:7){
#  file_name <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
# 
#  data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#  colnames(data)[1] <- "gene"
# 
#  # Sort the dataframe by avg_log2FC in descending order
#  data_sorted <- data[order(-data$avg_log2FC), ]
# 
#  # Select the top 20 genes from the first column
#  top_upregulated <- data_sorted[1:20, "gene"]
# 
#  # Sort the dataframe by avg_log2FC in descending order
#  data_sorted <- data[order(data$avg_log2FC), ]
# 
#  # Select the top 20 genes from the first column
#  top_downregulated <- data_sorted[1:20, "gene"]
# 
#  top_combined_genes <- c(top_upregulated, top_downregulated)
# 
#  # Create a heatmap for the cluster
# 
#  #Idents(object = lung.seurat.obj)  <- c('E0771_PBS', 'TF_PBS')
# 
#  Heatmap <- DoHeatmap(NDN_subset.seurat.obj, features = top_combined_genes, group.by = "group", assay = "RNA")
#  Heatmap <- Heatmap + theme(axis.text.y = element_text(size = 11, face = "bold"))
#  pdf_width <- 6
#  pdf_height <- 8
#  output_file2 <- file.path("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/", paste("Heatmap", i, ".pdf", sep = ""))
#  output_format <- "pdf"
#  ggsave(filename = output_file2, plot = Heatmap, device = output_format, width = pdf_width, height = pdf_height)
# 
# 
#  # Create a volcano plot for the cluster
# 
#  file_name2 <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i , "_NoFiltering.txt", sep = "")
# 
#  data2 <- read.table(file_name2, header = TRUE, sep = "\t")
# 
#  colnames(data2)[1] <- "gene"
# 
#  maxVal <- 1
#  minVal <- 1e-100
#  data2$p_val[data2$p_val < minVal] <- minVal
#  data2$p_val_adj[data2$p_val_adj < minVal] <- minVal
#  data2$p_val[is.na(data2$p_val)] <- maxVal
#  data2$p_val_adj[is.na(data2$p_val_adj)] <- maxVal
# 
#  data2 <- data2[order(data2$p_val), ]                                             #sort by p value
# 
#  data2$sig <- "NS"
#  data2$sig[data2$p_val <= 0.05] <- "p value <= 0.05"
#  data2$sig[data2$p_val_adj <= 0.05] <- "q value <= 0.05"
# 
#  maxCoord <- ceiling(max(data2$avg_log2FC))
#  minCoord <- floor(min(data2$avg_log2FC))
#  minPVal <- min(data2$p_val)
#  minCoord
#  maxCoord
#  maxYCoord <- ceiling(-log10(minPVal))
# 
#  topGenes <- head(data2, 20)                                                    # top 20 genes
#  topGenes
# 
#  genes = topGenes[0]
# 
#  title <-paste("C19_NDN_Markers_Cluster", i , sep = "")
# 
# 
#  volc = ggplot(data2, aes(avg_log2FC, -log10(p_val))) +
#    geom_point(aes(col=sig)) +                                                         # add points colored by significance
#    scale_color_manual(values=c("black", "pink", "red")) +
#    labs(title=title, x ="Log2FC", y = "-log10(p value)")  +
#    scale_x_continuous(limits = c(minCoord, maxCoord),breaks=pretty_breaks(n=11)) +                  # limits printing every other value
#    scale_y_continuous(limits = c(0, maxYCoord),breaks=pretty_breaks(n=11)) +                 # limits printing every value
# 
#    #     ###########################################################################################
# #   ## These next 3 lines will create bars at -1 and 1 with transparent shading between them ##
# #   ###########################################################################################
#  geom_vline(xintercept = 0, linetype="solid", color = "black", size=1.0) +
#    geom_vline(xintercept = 0,    linetype="solid", color = "black", size=1.0) +
#    annotate("rect", xmin=0, xmax=0, ymin=0, ymax=Inf, alpha=0.3, fill="black")  +
#    # Label the top 20 genes using geom_text_repel
#    geom_text_repel(data = topGenes, aes(label = gene), size = 3, box.padding = 0.5, max.overlaps = 50) +
# 
#    theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
#    #volc
#    output_file3 <- file.path("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/", paste("Volcano_Plot_Cluster", i , ".pdf", sep = ""))
#  pdf_width <- 6
#  pdf_height <- 8
#  output_format <- "pdf"
#  ggsave(filename = output_file3, plot = volc, device = output_format, width = pdf_width, height = pdf_height)
# }
# 
# 
# 
# #C19_NDN_VS_Healthy_NDN
# 
# for(i in 0:7){
# 
#   #THE OPTIONS ARE SET for no filtering of results.
#   seurat_object_sub<-subset(NDN_subset.seurat.obj, subset = seurat_clusters == i)
#   markers <- FindMarkers(seurat_object_sub, ident.1 = "C19_NDN",ident.2="Healthy_NDN", only.pos=FALSE, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, assay = "RNA")
#   fn<-paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster",i,"_NoFiltering.txt",sep="")
#   write.table(markers,fn,sep="\t", col.names = NA)
# 
#   #Filter markers by p_val_adj < 0.05
#   filtered_markers <- markers[markers$p_val_adj < 0.05, ]
# 
#   fn2 <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
#   write.table(filtered_markers, fn2, sep = "\t", col.names = NA)
# }
# 
# 
# for(i in 0:7){
#   file_name <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
# 
#   data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#   colnames(data)[1] <- "gene"
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(-data$avg_log2FC), ]
# 
#   # Select the top 20 genes from the first column
#   top_upregulated <- data_sorted[1:20, "gene"]
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(data$avg_log2FC), ]
# 
#   # Select the top 20 genes from the first column
#   top_downregulated <- data_sorted[1:20, "gene"]
# 
#   top_combined_genes <- c(top_upregulated, top_downregulated)
# 
#   # Create a heatmap for the cluster
# 
#   #Idents(object = lung.seurat.obj)  <- c('E0771_PBS', 'TF_PBS')
# 
#   Heatmap <- DoHeatmap(NDN_subset.seurat.obj, features = top_combined_genes, group.by = "group", assay = "RNA")
#   Heatmap <- Heatmap + theme(axis.text.y = element_text(size = 11, face = "bold"))
#   pdf_width <- 6
#   pdf_height <- 8
#   output_file2 <- file.path("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/", paste("Heatmap", i, ".pdf", sep = ""))
#   output_format <- "pdf"
#   ggsave(filename = output_file2, plot = Heatmap, device = output_format, width = pdf_width, height = pdf_height)
# 
# 
#   # Create a volcano plot for the cluster
# 
#   file_name2 <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i , "_NoFiltering.txt", sep = "")
# 
#   data2 <- read.table(file_name2, header = TRUE, sep = "\t")
# 
#   colnames(data2)[1] <- "gene"
# 
#   maxVal <- 1
#   minVal <- 1e-100
#   data2$p_val[data2$p_val < minVal] <- minVal
#   data2$p_val_adj[data2$p_val_adj < minVal] <- minVal
#   data2$p_val[is.na(data2$p_val)] <- maxVal
#   data2$p_val_adj[is.na(data2$p_val_adj)] <- maxVal
# 
#   data2 <- data2[order(data2$p_val), ]                                             #sort by p value
# 
#   data2$sig <- "NS"
#   data2$sig[data2$p_val <= 0.05] <- "p value <= 0.05"
#   data2$sig[data2$p_val_adj <= 0.05] <- "q value <= 0.05"
# 
#   maxCoord <- ceiling(max(data2$avg_log2FC))
#   minCoord <- floor(min(data2$avg_log2FC))
#   minPVal <- min(data2$p_val)
#   minCoord
#   maxCoord
#   maxYCoord <- ceiling(-log10(minPVal))
# 
#   topGenes <- head(data2, 20)                                                    # top 20 genes
#   topGenes
# 
#   genes = topGenes[0]
# 
#   title <-paste("C19_NDN_VS_Healthy_NDN_Cluster", i , sep = "")
# 
# 
#   volc = ggplot(data2, aes(avg_log2FC, -log10(p_val))) +
#     geom_point(aes(col=sig)) +                                                         # add points colored by significance
#     scale_color_manual(values=c("black", "pink", "red")) +
#     labs(title=title, x ="Log2FC", y = "-log10(p value)")  +
#     scale_x_continuous(limits = c(minCoord, maxCoord),breaks=pretty_breaks(n=11)) +                  # limits printing every other value
#     scale_y_continuous(limits = c(0, maxYCoord),breaks=pretty_breaks(n=11)) +                 # limits printing every value
# 
#     #     ###########################################################################################
#   #   ## These next 3 lines will create bars at -1 and 1 with transparent shading between them ##
#   #   ###########################################################################################
#   geom_vline(xintercept = 0, linetype="solid", color = "black", size=1.0) +
#     geom_vline(xintercept = 0,    linetype="solid", color = "black", size=1.0) +
#     annotate("rect", xmin=0, xmax=0, ymin=0, ymax=Inf, alpha=0.3, fill="black")  +
#     # Label the top 20 genes using geom_text_repel
#     geom_text_repel(data = topGenes, aes(label = gene), size = 3, box.padding = 0.5, max.overlaps = 50) +
# 
#     theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
#   #volc
#   output_file3 <- file.path("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/", paste("Volcano_Plot_Cluster", i , ".pdf", sep = ""))
#   pdf_width <- 6
#   pdf_height <- 8
#   output_format <- "pdf"
#   ggsave(filename = output_file3, plot = volc, device = output_format, width = pdf_width, height = pdf_height)
# }
# 
# 
# 
# #Non_C19_NDN_VS_Healthy_NDN
# 
# 
# 
# for(i in 0:7){
# 
#   #THE OPTIONS ARE SET for no filtering of results.
#   seurat_object_sub<-subset(NDN_subset.seurat.obj, subset = seurat_clusters == i)
#   markers <- FindMarkers(seurat_object_sub, ident.1 = "Non_C19_NDN",ident.2="Healthy_NDN", only.pos=FALSE, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, logfc.threshold = 0, assay = "RNA")
#   fn<-paste("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster",i,"_NoFiltering.txt",sep="")
#   write.table(markers,fn,sep="\t", col.names = NA)
# 
#   #Filter markers by p_val_adj < 0.05
#   filtered_markers <- markers[markers$p_val_adj < 0.05, ]
# 
#   fn2 <- paste("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
#   write.table(filtered_markers, fn2, sep = "\t", col.names = NA)
# }
# 
# 
# for(i in 0:7){
#   file_name <- paste("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i, "_Filtered.txt", sep = "")
# 
#   data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#   colnames(data)[1] <- "gene"
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(-data$avg_log2FC), ]
# 
#   # Select the top 20 genes from the first column
#   top_upregulated <- data_sorted[1:20, "gene"]
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(data$avg_log2FC), ]
# 
#   # Select the top 20 genes from the first column
#   top_downregulated <- data_sorted[1:20, "gene"]
# 
#   top_combined_genes <- c(top_upregulated, top_downregulated)
# 
#   # Create a heatmap for the cluster
# 
#   #Idents(object = lung.seurat.obj)  <- c('E0771_PBS', 'TF_PBS')
# 
#   Heatmap <- DoHeatmap(NDN_subset.seurat.obj, features = top_combined_genes, group.by = "group", assay = "RNA")
#   Heatmap <- Heatmap + theme(axis.text.y = element_text(size = 11, face = "bold"))
#   pdf_width <- 6
#   pdf_height <- 8
#   output_file2 <- file.path("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/", paste("Heatmap", i, ".pdf", sep = ""))
#   output_format <- "pdf"
#   ggsave(filename = output_file2, plot = Heatmap, device = output_format, width = pdf_width, height = pdf_height)
# 
# 
#   # Create a volcano plot for the cluster
# 
#   file_name2 <- paste("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i , "_NoFiltering.txt", sep = "")
# 
#   data2 <- read.table(file_name2, header = TRUE, sep = "\t")
# 
#   colnames(data2)[1] <- "gene"
# 
#   maxVal <- 1
#   minVal <- 1e-100
#   data2$p_val[data2$p_val < minVal] <- minVal
#   data2$p_val_adj[data2$p_val_adj < minVal] <- minVal
#   data2$p_val[is.na(data2$p_val)] <- maxVal
#   data2$p_val_adj[is.na(data2$p_val_adj)] <- maxVal
# 
#   data2 <- data2[order(data2$p_val), ]                                             #sort by p value
# 
#   data2$sig <- "NS"
#   data2$sig[data2$p_val <= 0.05] <- "p value <= 0.05"
#   data2$sig[data2$p_val_adj <= 0.05] <- "q value <= 0.05"
# 
#   maxCoord <- ceiling(max(data2$avg_log2FC))
#   minCoord <- floor(min(data2$avg_log2FC))
#   minPVal <- min(data2$p_val)
#   minCoord
#   maxCoord
#   maxYCoord <- ceiling(-log10(minPVal))
# 
#   topGenes <- head(data2, 20)                                                    # top 20 genes
#   topGenes
# 
#   genes = topGenes[0]
# 
#   title <-paste("Non_C19_NDN_VS_Healthy_NDN_Cluster", i , sep = "")
# 
# 
#   volc = ggplot(data2, aes(avg_log2FC, -log10(p_val))) +
#     geom_point(aes(col=sig)) +                                                         # add points colored by significance
#     scale_color_manual(values=c("black", "pink", "red")) +
#     labs(title=title, x ="Log2FC", y = "-log10(p value)")  +
#     scale_x_continuous(limits = c(minCoord, maxCoord),breaks=pretty_breaks(n=11)) +                  # limits printing every other value
#     scale_y_continuous(limits = c(0, maxYCoord),breaks=pretty_breaks(n=11)) +                 # limits printing every value
# 
#     #     ###########################################################################################
#   #   ## These next 3 lines will create bars at -1 and 1 with transparent shading between them ##
#   #   ###########################################################################################
#   geom_vline(xintercept = 0, linetype="solid", color = "black", size=1.0) +
#     geom_vline(xintercept = 0,    linetype="solid", color = "black", size=1.0) +
#     annotate("rect", xmin=0, xmax=0, ymin=0, ymax=Inf, alpha=0.3, fill="black")  +
#     # Label the top 20 genes using geom_text_repel
#     geom_text_repel(data = topGenes, aes(label = gene), size = 3, box.padding = 0.5, max.overlaps = 50) +
# 
#     theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
#   #volc
#   output_file3 <- file.path("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/", paste("Volcano_Plot_Cluster", i , ".pdf", sep = ""))
#   pdf_width <- 6
#   pdf_height <- 8
#   output_format <- "pdf"
#   ggsave(filename = output_file3, plot = volc, device = output_format, width = pdf_width, height = pdf_height)
# }
# 
# 
# 
# 
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #create table of DEGs counts at each comparison
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# # Create an empty data frame to store the results
# result_df <- data.frame(Round = numeric(0), All = numeric(0), Down = numeric(0), Up = numeric(0))
# 
# #C19_NDN_VS_Non_C19_NDN
# 
# for (i in 0:7){
# 
# 
#   file_name <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i , "_Filtered.txt", sep = "")
# 
#   data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#   colnames(data)[1] <- "gene"
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(-data$avg_log2FC), ]
# 
#   #count number of sig DEGs
#   All <- nrow(data_sorted)
# 
#   # Count the number of Down regulated DEGs
#   Down <- sum(data_sorted$avg_log2FC < 0)
# 
#   # Count the number of Up regulated DEGs
#   Up <- sum(data_sorted$avg_log2FC > 0)
# 
# 
#   # Append the results to the data frame
#   result_df <- rbind(result_df, data.frame(Round = i, All = All, Down = Down, Up = Up))
# 
# 
# }
# 
# # Specify the output file path
# output_file <- "./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Non_C19_NDN/DEG_count_table.txt"
# 
# # Write the data frame to a tab-delimited text file
# write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# ##########################################################################################################
# #C19_NDN_VS_Healthy_NDN
# 
# # Create an empty data frame to store the results
# result_df <- data.frame(Round = numeric(0), All = numeric(0), Down = numeric(0), Up = numeric(0))
# 
# 
# for (i in 0:7){
# 
# 
#  file_name <- paste("./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i , "_Filtered.txt", sep = "")
# 
#   data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#   colnames(data)[1] <- "gene"
# 
#  # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(-data$avg_log2FC), ]
# 
#   #count number of sig DEGs
#   All <- nrow(data_sorted)
# 
#   # Count the number of Down regulated DEGs
#   Down <- sum(data_sorted$avg_log2FC < 0)
# 
#   # Count the number of Up regulated DEGs
#   Up <- sum(data_sorted$avg_log2FC > 0)
# 
# 
#   # Append the results to the data frame
#   result_df <- rbind(result_df, data.frame(Round = i, All = All, Down = Down, Up = Up))
# 
# 
# }
# 
# # Specify the output file path
# output_file <- "./NDN_Analysis/my_rds/5DEGs/C19_NDN_VS_Healthy_NDN/DEG_count_table.txt"
# 
# # Write the data frame to a tab-delimited text file
# write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# ##########################################################################################################################
# # Create an empty data frame to store the results
# result_df <- data.frame(Round = numeric(0), All = numeric(0), Down = numeric(0), Up = numeric(0))
# 
# #Non_C19_NDN_VS_Healthy_NDN
# 
# for (i in 0:7){
# 
# 
#   file_name <- paste("./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i , "_Filtered.txt", sep = "")
# 
#   data <- read.table(file_name, header = TRUE, sep = "\t")
# 
#   colnames(data)[1] <- "gene"
# 
#   # Sort the dataframe by avg_log2FC in descending order
#   data_sorted <- data[order(-data$avg_log2FC), ]
# 
#   #count number of sig DEGs
#   All <- nrow(data_sorted)
# 
#   # Count the number of Down regulated DEGs
#   Down <- sum(data_sorted$avg_log2FC < 0)
# 
#   # Count the number of Up regulated DEGs
#   Up <- sum(data_sorted$avg_log2FC > 0)
# 
# 
#   # Append the results to the data frame
#   result_df <- rbind(result_df, data.frame(Round = i, All = All, Down = Down, Up = Up))
# 
# 
# }
# 
# # Specify the output file path
# output_file <- "./NDN_Analysis/my_rds/5DEGs/Non_C19_NDN_VS_Healthy_NDN/DEG_count_table.txt"
# 
# # Write the data frame to a tab-delimited text file
# write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pathway Enrichment in NDN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#C19_NDN_VS_Non_C19_NDN
# for (i in 0:7){
#   file_name <- paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames.txt", sep = "")
#   f<-read.table(file_name,sep="\t",header=TRUE)
#   df<-read.table(file_name,sep="\t",header=TRUE)
# 
#   pgostres <- gost(query = df$Gene,organism = "hsapiens", evcodes = TRUE)
#   #pgostres <- gost(query = df$EntrezID, organism = "hsapiens", evcodes = TRUE)
# 
# 
#  # The result is a named list where "result" is a data.frame with the enrichment analysis results
#  # and "meta" containing a named list with all the metadata for the query.
#   res <- as.data.frame(pgostres$result)
#   res$adjustedPValue<-p.adjust(res$p_value, method = "BH", n = length(res$p_value))
#   res <- res[order(res$source,res$adjustedPValue),]
#   res$parents <- vapply(res$parents, paste, collapse = ", ", character(1L))
#   res$evidence_codes <- vapply(res$evidence_codes, paste, collapse = ", ", character(1L))
#   res$intersection <- vapply(res$intersection, paste, collapse = ", ", character(1L))
#   fn<-paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.txt",sep="")
#   write.table(res,fn,row.names=FALSE,sep="\t")
# 
# 
#   df<-read.table(fn, sep="\t", header=TRUE)
#   df<-df[1:20,]
# 
#     #pdf(paste("./my_rds/6PathwayEnrichment/PT_VS_NT/PT_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""), width=10,height=10)
# 
#   #-adjustedPValue was used to plot the lowest adjP bars at the top(descending order)
#   plot = ggplot(df, aes(x=intersection_size, y=reorder(term_name, -adjustedPValue), fill=adjustedPValue)) +
#   geom_bar(stat="identity") +
#   xlab("Count")+
#   ylab("") +
#   guides(fill=guide_legend(title="Adjusted P Value")) +
#   scale_fill_gradientn(colours = c("red", "blue")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.title = element_text(size = 10),axis.title.y = element_text(size = 40))
# 
#   output_plot <- file.path(paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Non_C19_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""))
#   pdf_width <- 10
#   pdf_height <- 10
#   output_format <- "pdf"
#   ggsave(filename = output_plot, plot = plot, device = output_format, width = pdf_width, height = pdf_height)
# 
#     #dev.off()
# 
# }


#C19_NDN_VS_Healthy_NDN
# for (i in 0:7){
#   file_name <- paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames.txt", sep = "")
#   f<-read.table(file_name,sep="\t",header=TRUE)
#   df<-read.table(file_name,sep="\t",header=TRUE)
# 
#   pgostres <- gost(query = df$Gene,organism = "hsapiens", evcodes = TRUE)
#   #pgostres <- gost(query = df$EntrezID, organism = "hsapiens", evcodes = TRUE)
# 
# 
#   # The result is a named list where "result" is a data.frame with the enrichment analysis results
#   # and "meta" containing a named list with all the metadata for the query.
#   res <- as.data.frame(pgostres$result)
#   res$adjustedPValue<-p.adjust(res$p_value, method = "BH", n = length(res$p_value))
#   res <- res[order(res$source,res$adjustedPValue),]
#   res$parents <- vapply(res$parents, paste, collapse = ", ", character(1L))
#   res$evidence_codes <- vapply(res$evidence_codes, paste, collapse = ", ", character(1L))
#   res$intersection <- vapply(res$intersection, paste, collapse = ", ", character(1L))
#   fn<-paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.txt",sep="")
#   write.table(res,fn,row.names=FALSE,sep="\t")
# 
# 
#   df<-read.table(fn, sep="\t", header=TRUE)
#   df<-df[1:20,]
# 
#   #pdf(paste("./my_rds/6PathwayEnrichment/PT_VS_NT/PT_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""), width=10,height=10)
# 
#   #-adjustedPValue was used to plot the lowest adjP bars at the top(descending order)
#   plot = ggplot(df, aes(x=intersection_size, y=reorder(term_name, -adjustedPValue), fill=adjustedPValue)) +
#     geom_bar(stat="identity") +
#     xlab("Count")+
#     ylab("") +
#     guides(fill=guide_legend(title="Adjusted P Value")) +
#     scale_fill_gradientn(colours = c("red", "blue")) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           legend.title = element_text(size = 10),axis.title.y = element_text(size = 40))
# 
#   output_plot <- file.path(paste("./NDN_Analysis/my_rds/6PathwayEnrichment/C19_NDN_VS_Healthy_NDN/C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""))
#   pdf_width <- 10
#   pdf_height <- 10
#   output_format <- "pdf"
#   ggsave(filename = output_plot, plot = plot, device = output_format, width = pdf_width, height = pdf_height)
# 
#   #dev.off()
# 
# }
# 
# 
# #Non_C19_NDN_VS_Healthy_NDN
# 
# for (i in 0:7){
#   file_name <- paste("./NDN_Analysis/my_rds/6PathwayEnrichment/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames.txt", sep = "")
#   f<-read.table(file_name,sep="\t",header=TRUE)
#   df<-read.table(file_name,sep="\t",header=TRUE)
# 
#   pgostres <- gost(query = df$Gene,organism = "hsapiens", evcodes = TRUE)
#   #pgostres <- gost(query = df$EntrezID, organism = "hsapiens", evcodes = TRUE)
# 
# 
#   # The result is a named list where "result" is a data.frame with the enrichment analysis results
#   # and "meta" containing a named list with all the metadata for the query.
#   res <- as.data.frame(pgostres$result)
#   res$adjustedPValue<-p.adjust(res$p_value, method = "BH", n = length(res$p_value))
#   res <- res[order(res$source,res$adjustedPValue),]
#   res$parents <- vapply(res$parents, paste, collapse = ", ", character(1L))
#   res$evidence_codes <- vapply(res$evidence_codes, paste, collapse = ", ", character(1L))
#   res$intersection <- vapply(res$intersection, paste, collapse = ", ", character(1L))
#   fn<-paste("./NDN_Analysis/my_rds/6PathwayEnrichment/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.txt",sep="")
#   write.table(res,fn,row.names=FALSE,sep="\t")
# 
# 
#   df<-read.table(fn, sep="\t", header=TRUE)
#   df<-df[1:20,]
# 
#   #pdf(paste("./my_rds/6PathwayEnrichment/PT_VS_NT/PT_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""), width=10,height=10)
# 
#   #-adjustedPValue was used to plot the lowest adjP bars at the top(descending order)
#   plot = ggplot(df, aes(x=intersection_size, y=reorder(term_name, -adjustedPValue), fill=adjustedPValue)) +
#     geom_bar(stat="identity") +
#     xlab("Count")+
#     ylab("") +
#     guides(fill=guide_legend(title="Adjusted P Value")) +
#     scale_fill_gradientn(colours = c("red", "blue")) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           legend.title = element_text(size = 10),axis.title.y = element_text(size = 40))
# 
#   output_plot <- file.path(paste("./NDN_Analysis/my_rds/6PathwayEnrichment/Non_C19_NDN_VS_Healthy_NDN/Non_C19_NDN_Markers_Cluster", i, "_Filtered_GeneNames_gProfiler.pdf",sep=""))
#   pdf_width <- 10
#   pdf_height <- 10
#   output_format <- "pdf"
#   ggsave(filename = output_plot, plot = plot, device = output_format, width = pdf_width, height = pdf_height)
# 
#   #dev.off()
# 
# }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pseudotime in NDN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NDN_subset.seurat.obj <- readRDS("./NDN_Analysis/my_rds/covid.NDN.0.3.Markers.BTfilter.seurat.obj.rds")

NDN_subset.seurat.obj

NDN_sl <- slingshot(Embeddings(NDN_subset.seurat.obj, "umap"), clusterLabels = NDN_subset.seurat.obj$seurat_clusters,
                         start.clus = 4, stretch = 0)
pt <- slingPseudotime(NDN_sl)

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]

output_plot <- "./NDN_Analysis/my_rds/7Trajectory/NDN_Trajectory.pdf"
pdf(output_plot, width = 10, height = 7)
plot(reducedDims(NDN_sl), col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(NDN_sl), lwd = 2, col = 'black')
dev.off()






sessionInfo()






