library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dplyr)




covid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.seurat.obj.rds")

view(covid.seurat.obj)

unique(covid.seurat.obj@meta.data$sample)

covid.seurat.obj

covid.seurat.obj$MTPercent <- PercentageFeatureSet(covid.seurat.obj, pattern = '^MT-' )

view(covid.seurat.obj)

# If you want to create a new Seurat object with the modified data

covid.seurat.obj@meta.data$MTPercent[is.nan(covid.seurat.obj@meta.data$MTPercent)] <- 0

view(covid.seurat.obj)

violinplot1 <- VlnPlot(covid.seurat.obj, features = c("nCount_RNA", "nFeature_RNA", "MTPercent"), ncol = 3)
violinplot1

featurescatter <- FeatureScatter(covid.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE) + geom_smooth(method = 'lm')
featurescatter
covid.seurat.obj <- subset(covid.seurat.obj, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & MTPercent < 30 )



#I want to see how many cells left from a sample after filtering

covid.filtered.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.filtered.seurat.obj.rds")
covid.filtered.seurat.obj

sample_counts <- table(covid.filtered.seurat.obj$sample)
print(sample_counts)

#C19_P1_CD16_high      C19_P1_CD16_Int           C19_P1_NDN     C19_P2_CD16_high      C19_P2_CD16_Int           C19_P2_NDN 
#1408                 1942                 7570                 1142                 3687                13196 
#C19_P3_CD16_high      C19_P3_CD16_Int           C19_P3_NDN     C19_P4_CD16_high           C19_P4_NDN     C19_P5_CD16_high 
#12392                   75                 4262                 2379                 7856                  467 
#C19_P5_NDN     C19_P6_CD16_high           C19_P6_NDN     C19_P7_CD16_high      C19_P7_CD16_Int           C19_P7_NDN 
#3966                 4726                 3218                 3255                  581                 6369 
#Healthy_d1_NDN       Healthy_d2_NDN       Healthy_d3_NDN       Healthy_d4_NDN       Healthy_d5_NDN Non_C19_P1_CD16_high 
#11742                12822                  421                  924                  258                13015 
#Non_C19_P1_CD16_Int       Non_C19_P1_NDN Non_C19_P2_CD16_high       Non_C19_P2_NDN Post_C19_1_CD16_high       Post_C19_1_NDN 
#11543                10960                  328                 3078                  109                 7398 

view(covid.filtered.seurat.obj)

#remove a sample from the object C19_P3_CD16_Int before integration. Because this sample has <100 cells

# Define the sample name you want to remove
sample_to_remove <- c("C19_P3_CD16_Int")  # Replace "SampleName" with the name of the sample you want to remove

# Find cells belonging to the specified sample
cells_to_remove <- which(covid.filtered.seurat.obj@meta.data$sample == sample_to_remove)

# Create a new Seurat object excluding the specified sample
covid.filtered.seurat.obj.Subset <- subset(covid.filtered.seurat.obj, cells = -cells_to_remove)

covid.filtered.seurat.obj
#An object of class Seurat 
#37143 features across 151089 samples within 1 assay 
#Active assay: RNA (37143 features, 2000 variable features)

covid.filtered.seurat.obj.Subset
#An object of class Seurat 
#37143 features across 151014 samples within 1 assay 
#Active assay: RNA (37143 features, 2000 variable features)



covid.filtered.integrated.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.filtered.integrated.seurat.obj.rds")
covid.filtered.integrated.seurat.obj


covid.filtered.integrated.seurat.obj <- RunUMAP(object = covid.filtered.integrated.seurat.obj, dims = 1:30)


p3 <- DimPlot(covid.filtered.integrated.seurat.obj, reduction = 'umap', group.by = 'sample')
p4 <- DimPlot(covid.filtered.integrated.seurat.obj, reduction = 'umap', group.by = 'group', cols = c('red','green','blue','yellow','violet','cyan','orange','gray','pink'))
p5 <- DimPlot(covid.filtered.integrated.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)

covid.filter3.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.filtered.seurat.obj.EmptyDrops.rds")
covid.filter3.seurat.obj
#An object of class Seurat 
#37143 features across 120981 samples within 1 assay 
#Active assay: RNA (37143 features, 2000 variable features)
view(covid.filter3.seurat.obj)

covid.filter3.harmony.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA//Harmony/covid.filtered.seurat.HarmonyIntegration.obj.EmptyDrops.rds")
covid.filter3.harmony.seurat.obj
#An object of class Seurat 
#37143 features across 120981 samples within 1 assay 
#Active assay: RNA (37143 features, 2000 variable features)
#3 dimensional reductions calculated: pca, harmony, umap
view(covid.filter3.harmony.seurat.obj)



library(Seurat)

#dotplot of genes of interest

covid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.Allmarkers.seurat.BTfilter.obj.rds")

view(covid.seurat.obj)

#Idents(object = NDNcovid.seurat.obj) <- "RNA_snn_res_0.3"
Idents(object = covid.seurat.obj) <- "seurat_clusters"
DefaultAssay(covid.seurat.obj) <- "RNA"


# List of top 10 genes
top3_genes <- c("CYTH4", "SYAP1", "PGGHG", "IFIT2", "IFIT1", "IFIT3", "LTF", "LCN2", "CAMP", "AZU1", "DEFA3", "MPO", "DEPRECATED-ENSG00000276085", "CCL3", "CCL4", "ID2", "DEPRECATED-ENSG00000115738", "PTGS2", "LY6E", "ISG15", "RSAD2", "PRF1", "CLC", "GNLY", "NRGN", "PPBP", "ITGA2B", "DEFA3", "DEFA4", "LYZ", "CAMP", "MMP8", "CRISP3", "VCAN", "THBS1", "MAFB") 

# Ensure that the 'seurat_clusters' column is a factor
covid.seurat.obj$seurat_clusters <- factor(covid.seurat.obj$seurat_clusters)

any(is.na(covid.seurat.obj$seurat_clusters))


# Create dot plots for each top gene
dot_plot <- DotPlot(object = covid.seurat.obj, features = top3_genes,
                    group.by = "seurat_clusters",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
 ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/All_Res0.3/Top3Markers_dotplot_seuratcluster.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)

# Create dot plots for each top gene
dot_plot <- DotPlot(object = covid.seurat.obj, features = top3_genes,
                    group.by = "group",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
 ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/All_Res0.3/Top3Markers_dotplot_group.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)


#Bycelltype
#NDN
NDNcovid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.NDN.0.3.Markers.seurat.obj.rds")
view(NDNcovid.seurat.obj)

#Idents(object = NDNcovid.seurat.obj) <- "RNA_snn_res_0.3"
Idents(object = NDNcovid.seurat.obj) <- "seurat_clusters"
DefaultAssay(NDNcovid.seurat.obj) <- "RNA"


# List of top 10 genes
top3_genes <- c("MMP9", "CR1", "CKAP4", "RNF149", "STK4", "AMPD2", "RSAD2", "GBP1", "ISG15", "PLIN5", "PLIN4", "ITGAX", "CCL4", "CCL3", "CCL4L2", "MT-CO3", "MT-ND4", "MT-CO1", "LTF", "LCN2", "ARG1", "AZU1", "DEFA3", "MPO", "CST3", "CD68", "CD14", "PRF1", "GNLY", "JCHAIN") 

# Ensure that the 'seurat_clusters' column is a factor
NDNcovid.seurat.obj$seurat_clusters <- factor(NDNcovid.seurat.obj$seurat_clusters)



# Create dot plots for each top gene
dot_plot <- DotPlot(object = NDNcovid.seurat.obj, features = top3_genes,
                    group.by = "seurat_clusters",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                    

pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/NDN/Top3Markers_dotplot_seuratcluster.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)

# Create dot plots for each top gene
dot_plot <- DotPlot(object = NDNcovid.seurat.obj, features = top3_genes,
                    group.by = "group",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/NDN/Top3Markers_dotplot_group.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)


#CD16_high
CD16_highcovid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.CD16_high.0.3.Markers.seurat.obj.rds")
#view(CD16_highcovid.seurat.obj)

#Idents(object = NDNcovid.seurat.obj) <- "RNA_snn_res_0.3"
Idents(object = CD16_highcovid.seurat.obj) <- "seurat_clusters"
DefaultAssay(CD16_highcovid.seurat.obj) <- "RNA"


# List of top 10 genes
top3_genes <- c("MMP9", "S100A12", "S100A8", "IFI44", "IFI44L", "NAPA", "ADGRE2", "PPIF", "AC099489.1", "DEPRECATED-ENSG00000276085", "CCL4L2", "CCL3L1", "NRGN", "PPBP", "ITGB3", "LTF", "LCN2", "CAMP", "PRF1", "NKG7", "GNLY", "VCAN", "CCL3", "THBS1") 

# Ensure that the 'seurat_clusters' column is a factor
CD16_highcovid.seurat.obj$seurat_clusters <- factor(CD16_highcovid.seurat.obj$seurat_clusters)



# Create dot plots for each top gene
dot_plot <- DotPlot(object = CD16_highcovid.seurat.obj, features = top3_genes,
                    group.by = "seurat_clusters",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/CD16_high/Top3Markers_dotplot_seuratcluster.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)

# Create dot plots for each top gene
dot_plot <- DotPlot(object = CD16_highcovid.seurat.obj, features = top3_genes,
                    group.by = "group",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/CD16_high/Top3Markers_dotplot_group.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)



#CD16_int
CD16_intcovid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/covid.CD16_int.0.3.Markers.seurat.obj.rds")
#view(CD16_highcovid.seurat.obj)

#Idents(object = NDNcovid.seurat.obj) <- "RNA_snn_res_0.3"
Idents(object = CD16_intcovid.seurat.obj) <- "seurat_clusters"
DefaultAssay(CD16_intcovid.seurat.obj) <- "RNA"


# List of top 10 genes
top3_genes <- c("PLBD1", "CRISP3", "LTF", "CAMP", "CHI3L1", "ABCA13", "MMP9", "PTGS2", "FPR1", "SMPD3", "ADGRG5", "LFNG", "AZU1", "MPO", "ELANE", "HIST1H1B", "RRM2", "NUSAP1", "CLEC7A", "NAMPT", "EGR3", "VCAN", "AHNAK", "PRF1") 

# Ensure that the 'seurat_clusters' column is a factor
CD16_intcovid.seurat.obj$seurat_clusters <- factor(CD16_intcovid.seurat.obj$seurat_clusters)



# Create dot plots for each top gene
dot_plot <- DotPlot(object = CD16_intcovid.seurat.obj, features = top3_genes,
                    group.by = "seurat_clusters",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/CD16_int/Top3Markers_dotplot_seuratcluster.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)

# Create dot plots for each top gene
dot_plot <- DotPlot(object = CD16_intcovid.seurat.obj, features = top3_genes,
                    group.by = "group",  # Replace with the actual group column name
                    assay = "RNA",
                    #cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/my_rds/4Markers/CD16_int/Top3Markers_dotplot_group.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)


#################################################
#stacked bar plot
#################################################

#NDN

data <- read.csv("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_clusters_percent.csv")
df <- data.frame(data)
df

Cluster_ID <- c()
for (id in df$Cluster_ID)
  #replicate the Cluste_ID 3 times becaus ewe have 3 subcategoreies in each cluster
  Cluster_ID <- c(Cluster_ID, rep(id, 3))

Cluster_ID

Cells <- colnames(df[, 2:4])
Cells

n_replications <- nrow(df)
Cells <- rep(Cells, n_replications)
Cells

# get the columns that contains the opinion values
values <- df[2:4]
values

values <- t(values) 
values

values <- as.vector(values)
values


final_df = data.frame(Cluster_ID, Cells, values)
colnames(final_df)[3]  <- "Cell Percent"
head(final_df)


library(ggplot2)

# Ensure Cluster_ID is a factor with levels 0 to 9
final_df$Cluster_ID <- factor(final_df$Cluster_ID, levels = 0:9)

stackedbarplot <- ggplot(final_df, aes(fill=Cells, y=`Cell Percent`, x=Cluster_ID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
        scale_x_discrete(limits = as.character(0:9)) 

ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_stackedbarplot.pdf", plot = stackedbarplot, width = 8, height = 6, units = "in")

# ggplot(final_df, aes(x=Cluster_ID, y=`Cell Percent`, fill=cells)) + 
#   geom_col(position="dodge") + 
#   theme(axis.text.x = element_text(angle = 45, margin = margin(t=30, "pt")))

barplot <- ggplot(final_df, aes(x=Cluster_ID, y=`Cell Percent`, fill=cells)) + 
  geom_col(position="dodge") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(0:9)) 


ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_barplot.pdf", plot = barplot, width = 8, height = 6, units = "in")


#high & int

data <- read.csv("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/CD16_int_clusters_percent.csv")
df <- data.frame(data)
df

Cluster_ID <- c()
for (id in df$Cluster_ID)
  #replicate the Cluste_ID 3 times becaus ewe have 3 subcategoreies in each cluster
  Cluster_ID <- c(Cluster_ID, rep(id, 2))

Cluster_ID

Cells <- colnames(df[, 2:3])
Cells

n_replications <- nrow(df)
Cells <- rep(Cells, n_replications)
Cells

# get the columns that contains the opinion values
values <- df[2:3]
values

values <- t(values) 
values

values <- as.vector(values)
values


final_df = data.frame(Cluster_ID, Cells, values)
colnames(final_df)[3]  <- "Cell Percent"
head(final_df)


library(ggplot2)

# Ensure Cluster_ID is a factor with levels 0 to 7
final_df$Cluster_ID <- factor(final_df$Cluster_ID, levels = 0:7)

stackedbarplot <- ggplot(final_df, aes(fill=Cells, y=`Cell Percent`, x=Cluster_ID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(0:7)) 

ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/CD16_int_stackedbarplot.pdf", plot = stackedbarplot, width = 8, height = 6, units = "in")

# ggplot(final_df, aes(x=Cluster_ID, y=`Cell Percent`, fill=cells)) + 
#   geom_col(position="dodge") + 
#   theme(axis.text.x = element_text(angle = 45, margin = margin(t=30, "pt")))

barplot <- ggplot(final_df, aes(x=Cluster_ID, y=`Cell Percent`, fill=Cells)) + 
  geom_col(position="dodge") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(0:7)) 


ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/CD16_int_barplot.pdf", plot = barplot, width = 8, height = 6, units = "in")

#######################################################################
#pie chart
#######################################################################

#Test1
# Load necessary library
library(ggplot2)

# Create the data frame
df <- data.frame(
  `Cluster_ID` = c(0, 1, 2, 3, 4, 5, 6, 7),
  `C19_NDN` = c(19.17942, 35.35656, 18.75611, 14.99512, 6.222729, 1.794204, 2.61804, 1.077825),
  `Healthy_NDN` = c(60.70331, 13.21746, 14.519, 2.332256, 0.521423, 4.919159, 1.204527, 2.582862),
  `Non_C19_NDN` = c(7.91165, 46.18448, 14.49262, 26.47778, 0.470716, 1.167738, 2.136327, 1.158686)
)

bp<- ggplot(df, aes(x="", y=C19_NDN, fill=Cluster_ID))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + coord_polar("y", start=0) + geom_text(aes(label = C19_NDN, size=5))
pie

########################correct code######################
# Load the necessary library
library(ggplot2)
library(ggrepel)
library(dplyr)

# Create the data frame
df <- data.frame(
  Cluster_ID = c(0, 1, 2, 3, 4, 5, 6, 7),
  C19_NDN = c(19.17942, 35.35656, 18.75611, 14.99512, 6.222729, 1.794204, 2.61804, 1.077825),
  Healthy_NDN = c(60.70331, 13.21746, 14.519, 2.332256, 0.521423, 4.919159, 1.204527, 2.582862),
  Non_C19_NDN = c(7.91165, 46.18448, 14.49262, 26.47778, 0.470716, 1.167738, 2.136327, 1.158686)
)

# Define pastel colors
pastel_colors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC")

###### Create the pie chart for C19_NDN

# Calculate the position of labels
df <- df %>% 
  mutate(csum = rev(cumsum(rev(C19_NDN))), 
         pos = C19_NDN/2 + lead(csum, 1),
         pos = if_else(is.na(pos), C19_NDN/2, pos))

bp_c19 <- ggplot(df, aes(x="", y=C19_NDN, fill=factor(Cluster_ID))) +
  geom_bar(width=1, stat="identity", color="white") +
  coord_polar("y", start=0) +
  theme_void() +  # Remove background, grid, and axis marks
  scale_fill_manual(values=pastel_colors) +  # Use pastel colors
  geom_text_repel(aes(y = pos, label = paste0(round(C19_NDN, 3), "%")), size = 5, fontface = "bold", nudge_x = 0.7, segment.color = NA, show.legend = FALSE) +  # Add labels without overlapping
  #geom_label_repel(aes(y = pos, label = paste0(round(C19_NDN, 3), "%")), size = 4.5, nudge_x = 1, show.legend = FALSE) +
  labs(fill="Cluster ID", title="C19 NDN(%) across clusters", size = 12, face = "bold")+  # Add legend and title
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 12, face = "bold"))

# Print the pie chart
print(bp_c19)

#ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/C19_NDN_pie.pdf", plot = bp_c19) 

pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/C19_NDN_pie.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = bp_c19, device = output_format, width = pdf_width, height = pdf_height)


##### Create the pie chart for Healthy_NDN

# Calculate the position of labels
df <- df %>% 
  mutate(csum = rev(cumsum(rev(Healthy_NDN))), 
         pos = Healthy_NDN/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Healthy_NDN/2, pos))

bp_Healthy <- ggplot(df, aes(x="", y=Healthy_NDN, fill=factor(Cluster_ID))) +
  geom_bar(width=1, stat="identity", color="white") +
  coord_polar("y", start=0) +
  theme_void() +  # Remove background, grid, and axis marks
  scale_fill_manual(values=pastel_colors) +  # Use pastel colors
  geom_text_repel(aes(y = pos, label = paste0(round(Healthy_NDN, 3), "%")), size = 5, fontface = "bold", nudge_x = 0.7, segment.color = NA, show.legend = FALSE) +  # Add labels without overlapping
  #geom_label_repel(aes(y = pos, label = paste0(round(C19_NDN, 3), "%")), size = 4.5, nudge_x = 1, show.legend = FALSE) +
  labs(fill="Cluster ID", title="Healthy NDN(%) across clusters", size = 12, face = "bold")+  # Add legend and title
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 12, face = "bold"))

# Print the pie chart
print(bp_Healthy)

#ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/C19_NDN_pie.pdf", plot = bp_c19) 

pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Healthy_NDN_pie.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = bp_Healthy, device = output_format, width = pdf_width, height = pdf_height)


##### Create the pie chart for Non_C19_NDN

# Calculate the position of labels
df <- df %>% 
  mutate(csum = rev(cumsum(rev(Non_C19_NDN))), 
         pos = Non_C19_NDN/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Non_C19_NDN/2, pos))

bp_nonc19 <- ggplot(df, aes(x="", y=Non_C19_NDN, fill=factor(Cluster_ID))) +
  geom_bar(width=1, stat="identity", color="white") +
  coord_polar("y", start=0) +
  theme_void() +  # Remove background, grid, and axis marks
  scale_fill_manual(values=pastel_colors) +  # Use pastel colors
  geom_text_repel(aes(y = pos, label = paste0(round(Non_C19_NDN, 3), "%")), size = 5, fontface = "bold", nudge_x = 0.7, segment.color = NA, show.legend = FALSE) +  # Add labels without overlapping
  #geom_label_repel(aes(y = pos, label = paste0(round(C19_NDN, 3), "%")), size = 4.5, nudge_x = 1, show.legend = FALSE) +
  labs(fill="Cluster ID", title="Non_C19_NDN(%) across clusters", size = 12, face = "bold")+  # Add legend and title
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 12, face = "bold"))

# Print the pie chart
print(bp_nonc19)

#ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/C19_NDN_pie.pdf", plot = bp_c19) 

pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Non_C19_NDN_pie.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = bp_nonc19, device = output_format, width = pdf_width, height = pdf_height)

###################### Simple Pie Chart
slices <- c(19.17942, 35.35656, 18.75611, 14.99512, 6.222729, 1.794204, 2.61804, 1.077825)
labels <- c("0", "1", "2", "3", "4", "5", "6", "7")
C19_NDN_pie <- pie(slices, labels, main="C19 NDN distribution across clusters", col=pastel_colors, border="white")

  



#################################################
#stacked bar plot
#################################################

#C19_NDN

df <- data.frame(
  C19_Patient_ID = c(1, 2, 3, 4, 5, 6, 7),
  Cluster0 = c(13.2108,	16.60955,	7.22777,	23.51773,	35.52732,	0,	0),
  Cluster1 = c(13.9575,	29.13716,	30.26629,	47.96188,	46.91233,	57.77778,	20.83333),
  Cluster2 = c(0.804136,	22.83147,	42.55825,	9.502382,	5.006353,	0,	44.79167),
  Cluster3 = c(8.328547,	22.81624,	10.17594,	12.2684,	2.592122,	0,	7.291667),
  Cluster4 = c(58.4147,	4.226639,	5.991441,	0.913182,	0.355781,	6.666667,	1.041667),
  Cluster5 = c(2.125215,	0.769172,	0.42796,	2.567496,	4.98094,	6.666667,	2.083333),
  Cluster6 = c(1.665709,	2.520752,	2.78174,	2.170461,	3.939009,	6.666667,	5.208333),
  Cluster7 = c(1.493395,	1.089026,	0.570613,	1.098465,	0.68615,	22.22222,	18.75)
)

C19_Patient_ID <- c()
for (id in df$C19_Patient_ID)
  #replicate the Cluste_ID 3 times becaus ewe have 3 subcategoreies in each cluster
  C19_Patient_ID <- c(C19_Patient_ID, rep(id, 8))

C19_Patient_ID

Clusters <- colnames(df[, 2:9])
Clusters

n_replications <- nrow(df)
Clusters <- rep(Clusters, n_replications)
Clusters

# get the columns that contains the opinion values
values <- df[2:9]
values

values <- t(values) 
values

values <- as.vector(values)
values


final_df = data.frame(C19_Patient_ID, Clusters, values)
colnames(final_df)[3]  <- "NDN Percent"
head(final_df)
final_df

library(ggplot2)

# Ensure Cluster_ID is a factor with levels 0 to 9
final_df$C19_Patient_ID <- factor(final_df$C19_Patient_ID, levels = 1:7)

stackedbarplot <- ggplot(final_df, aes(fill=Clusters, y=`NDN Percent`, x=C19_Patient_ID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(1:7)) 

ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/C19_patient_stackedbarplot.pdf", plot = stackedbarplot, width = 8, height = 6, units = "in")


#Healthy_NDN

df <- data.frame(
  Healthy_Patient_ID = c(1, 2, 3, 4, 5),
  Cluster0 = c(66.90819,	54.86035,	72.72727,	72.72727,	69.49153),
  Cluster1 = c(15.32691,	11.46045,	3.636364,	0,	1.694915),
  Cluster2 = c(4.628761,	23.67764,	3.636364,	5.194805,	13.55932),
  Cluster3 = c(3.350098,	1.435481,	0,	0,	0),
  Cluster4 = c(0.639332,	0.413481,	0,	1.298701,	0),
  Cluster5 = c(4.679908,	5.086597,	9.090909,	10.38961,	5.084746),
  Cluster6 = c(0.699003,	1.661726,	3.636364,	1.298701,	0),
  Cluster7 = c(3.767795,	1.404275,	7.272727,	9.090909,	10.16949)
)

Healthy_Patient_ID <- c()
for (id in df$Healthy_Patient_ID)
  #replicate the Cluste_ID 3 times becaus ewe have 3 subcategoreies in each cluster
  Healthy_Patient_ID <- c(Healthy_Patient_ID, rep(id, 8))

Healthy_Patient_ID

Clusters <- colnames(df[, 2:9])
Clusters

n_replications <- nrow(df)
Clusters <- rep(Clusters, n_replications)
Clusters

# get the columns that contains the opinion values
values <- df[2:9]
values

values <- t(values) 
values

values <- as.vector(values)
values


final_df = data.frame(Healthy_Patient_ID, Clusters, values)
colnames(final_df)[3]  <- "NDN Percent"
head(final_df)
final_df

library(ggplot2)

# Ensure Cluster_ID is a factor with levels 0 to 9
final_df$Healthy_Patient_ID <- factor(final_df$Healthy_Patient_ID, levels = 1:5)

stackedbarplot <- ggplot(final_df, aes(fill=Clusters, y=`NDN Percent`, x=Healthy_Patient_ID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(1:5)) 

ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Healthy_patient_stackedbarplot.pdf", plot = stackedbarplot, width = 8, height = 6, units = "in")



#Non_C19_NDN

df <- data.frame(
  Non_C19_Patient_ID = c(1, 2),
  Cluster0 = c(7.863014,	13.40206),
  Cluster1 = c(46.05479,	60.82474),
  Cluster2 = c(14.59361,	3.092784),
  Cluster3 = c(26.6758,	4.123711),
  Cluster4 = c(0.474886,	0),
  Cluster5 = c(1.09589,	9.278351),
  Cluster6 = c(2.091324,	7.216495),
  Cluster7 = c(1.150685,	2.061856)
)

Non_C19_Patient_ID <- c()
for (id in df$Non_C19_Patient_ID)
  #replicate the Cluste_ID 3 times becaus ewe have 3 subcategoreies in each cluster
  Non_C19_Patient_ID <- c(Non_C19_Patient_ID, rep(id, 8))

Non_C19_Patient_ID

Clusters <- colnames(df[, 2:9])
Clusters

n_replications <- nrow(df)
Clusters <- rep(Clusters, n_replications)
Clusters

# get the columns that contains the opinion values
values <- df[2:9]
values

values <- t(values) 
values

values <- as.vector(values)
values


final_df = data.frame(Non_C19_Patient_ID, Clusters, values)
colnames(final_df)[3]  <- "NDN Percent"
head(final_df)
final_df

library(ggplot2)

# Ensure Cluster_ID is a factor with levels 0 to 9
final_df$Non_C19_Patient_ID <- factor(final_df$Non_C19_Patient_ID, levels = 1:2)

stackedbarplot <- ggplot(final_df, aes(fill=Clusters, y=`NDN Percent`, x=Non_C19_Patient_ID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_x_discrete(limits = as.character(1:2)) 

ggsave("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Non_C19_patient_stackedbarplot.pdf", plot = stackedbarplot, width = 8, height = 6, units = "in")

###########################################################
#dotplot
###########################################################
library(Seurat)

#NDN
NDNcovid.seurat.obj <- readRDS("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/covid.NDN.0.3.Markers.BTfilter.seurat.obj.rds")
NDNcovid.seurat.obj

#Idents(object = NDNcovid.seurat.obj) <- "RNA_snn_res_0.3"
Idents(object = NDNcovid.seurat.obj) <- "seurat_clusters"
DefaultAssay(NDNcovid.seurat.obj) <- "RNA"

# Ensure that the 'seurat_clusters' column is a factor with unique levels
#NDNcovid.seurat.obj$seurat_clusters <- factor(NDNcovid.seurat.obj$seurat_clusters, levels = unique(NDNcovid.seurat.obj$seurat_clusters))

# List of top 10 genes
top3_genes <- c("AMPD2", "DOCK5", "ATG2A", "IL1B", "TNFAIP3", "IL1RN", "GBP5", "GBP1", "RSAD2", "MMP9", "CD177", "S100A12", "AZU1", "LCN2", "DEFA3", "MT-CO3", "MT-ND4", "MT-CO1", "ASTL", "DUSP2", "TNFAIP3", "CST3", "CD68", "MAFB") 



# Create dot plots for each top gene
dot_plot <- DotPlot(object = NDNcovid.seurat.obj, features = unique(top3_genes),
                    group.by = "seurat_clusters",  # Replace with the actual group column name
                    assay = "RNA",
                    cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Top3Markers_dotplot_seuratcluster.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)

# Create dot plots for each top gene
dot_plot <- DotPlot(object = NDNcovid.seurat.obj, features = unique(top3_genes),
                    group.by = "group",  # Replace with the actual group column name
                    assay = "RNA",
                    cols = c("blue", "red"),  # Customize the colors for the groups
                    #min.cutoff = 0,
                    #max.cutoff = 4,
                    scale.by = "size"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf_width <- 10
pdf_height <- 7
output_file2 <- file.path("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA/NDN_Analysis/Top3Markers_dotplot_group.pdf")
output_format <- "pdf"
ggsave(filename = output_file2, plot = dot_plot, device = output_format, width = pdf_width, height = pdf_height)









library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)  # For enframe function

# List of top 3 genes for each cluster, allowing duplicates
top3_genes <- c("AMPD2", "DOCK5", "ATG2A", "IL1B", "TNFAIP3", "IL1RN", "GBP5", "GBP1", "RSAD2", "MMP9", "CD177", "S100A12", "AZU1", "LCN2", "DEFA3", "MT-CO3", "MT-ND4", "MT-CO1", "ASTL", "DUSP2", "TNFAIP3", "CST3", "CD68", "MAFB")

# Ensure the 'seurat_clusters' column is a factor
NDNcovid.seurat.obj$seurat_clusters <- factor(NDNcovid.seurat.obj$seurat_clusters)

# Extract data for unique genes
dot_data <- DotPlot(object = NDNcovid.seurat.obj, features = unique(top3_genes), group.by = "seurat_clusters", assay = "RNA", scale.by = "size")$data

# Replicate the data for duplicated genes
replicated_dot_data <- top3_genes %>% 
  enframe(name = NULL, value = "features.plot") %>%
  left_join(dot_data, by = "features.plot")

# Plot using ggplot2
dot_plot <- ggplot(replicated_dot_data, aes(x = id, y = features.plot, size = avg.exp.scaled, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "#B3CDE3", high = "#FBB4AE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  guides(size = guide_legend(title = "Scaled Expression"), color = guide_colorbar(title = "Scaled Expression"))

# Print the dot plot
print(dot_plot)
