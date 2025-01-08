
#March 2024 -Kalpani de Silva

library(Seurat)
#library(SeuratData)
library(ggplot2)
#library(patchwork)
library(dplyr)
library(Matrix)
library(harmony)
library(gridExtra)


#setwd("C:/Kalpani/KY-INBRE PostDoc/Single_Cell/Jun Yan/Fixed_RNA")
setwd("/bio/home/kkdesi01/KBRIN0XXX-Yan-FixedRNA")

#################
covid.filtered.seurat.obj <- readRDS("./NDN_Analysis/my_rds/covid.filtered.Bfilter.seurat.obj.EmptyDrops.rds")

##############################################################
#normalize, harmonize, cluster
covid.filtered.seurat.obj <- NormalizeData(covid.filtered.seurat.obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
covid.filtered.seurat.obj <- RunHarmony(covid.filtered.seurat.obj, group.by.vars = "sample", plot_convergence = TRUE, assay.use = "RNA", project.dim = F)
covid.filtered.seurat.obj <- RunUMAP(covid.filtered.seurat.obj, reduction = "harmony", dims = 1:30)

##############################################################0.5
covid.filtered.seurat.obj <- FindNeighbors(covid.filtered.seurat.obj, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.5)

#visualize
pdf("./NDN_Analysis/my_rds/2Integration/seuratclusters_EmptyDrops_BTfilter_Res0.5.pdf",width=30,height=20)
DimPlot(covid.filtered.seurat.obj, group.by = c("sample","group","celltype","seurat_clusters"), ncol = 2)
dev.off()
#13 clusters

##############################################################0.4
covid.filtered.seurat.obj <- FindNeighbors(covid.filtered.seurat.obj, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.4)

#visualize
pdf("./NDN_Analysis/my_rds/2Integration/seuratclusters_EmptyDrops_BTfilter_Res0.4.pdf",width=30,height=20)
DimPlot(covid.filtered.seurat.obj, group.by = c("sample","group","celltype","seurat_clusters"), ncol = 2)
dev.off()
#12 clusters

##############################################################0.3
covid.filtered.seurat.obj <- FindNeighbors(covid.filtered.seurat.obj, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.3)

#visualize
pdf("./NDN_Analysis/my_rds/2Integration/seuratclusters_EmptyDrops_BTfilter_Res0.3.pdf",width=30,height=20)
DimPlot(covid.filtered.seurat.obj, group.by = c("sample","group","celltype","seurat_clusters"), ncol = 2)
dev.off()
#8 clusters



p3 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'sample')
p31 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'group')
p32 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'celltype')
p33 <- DimPlot(covid.filtered.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, label.size = 6)
grid3 <- grid.arrange(p3, p31, p32, p33, ncol = 2, nrow = 2)

# Save the plot to the specified file
output_file6 <- "./NDN_Analysis/my_rds/2Integration/seuratclusters_EmptyDrops_BTfilter_Res0.3_withLabels.pdf"
output_format <- "pdf"
pdf_width <- 30
pdf_height <- 25
ggsave(filename = output_file6, plot = grid3, device = output_format, width = pdf_width, height = pdf_height) 

#pdf of only seurat clusters
pdf("./NDN_Analysis/my_rds/2Integration/Onlyseuratclusters_EmptyDrops_BTfilter_Res0.3.pdf",width=10,height=7)
DimPlot(covid.filtered.seurat.obj, group.by = "seurat_clusters", label = TRUE, label.size = 6)
dev.off()

#how many cells
sample_counts <- table(covid.filtered.seurat.obj$sample)
print(sample_counts)
#C19_P1_NDN     C19_P2_NDN     C19_P3_NDN     C19_P4_NDN     C19_P5_NDN 
#1741          13131           4206           7556           3935 
#C19_P6_NDN     C19_P7_NDN Healthy_d1_NDN Healthy_d2_NDN Healthy_d3_NDN 
#45             96          11731          12818             55 
#Healthy_d4_NDN Healthy_d5_NDN Non_C19_P1_NDN Non_C19_P2_NDN 
#77             59          10950             97                

cluster_count <- table(covid.filtered.seurat.obj$RNA_snn_res.0.3)
print(cluster_count)

#0     1     2     3     4     5     6     7 
#21782 19230 10953  8107  2092  1897  1338  1098   


# Identify clusters and groups
clusters <- levels(covid.filtered.seurat.obj$RNA_snn_res.0.3)
groups <- unique(covid.filtered.seurat.obj$group)

# Initialize a matrix to store the cell counts
cell_counts <- matrix(0, nrow = length(clusters), ncol = length(groups),
                      dimnames = list(clusters, groups))

# Count the cells in each cluster per group
for (cluster in clusters) {
  for (group in groups) {
    cell_counts[cluster, group] <- sum(covid.filtered.seurat.obj$RNA_snn_res.0.3 == cluster &
                                         covid.filtered.seurat.obj$group == group)
  }
}

print(cell_counts)

#C19_NDN Healthy_NDN Non_C19_NDN
#0    5890       15018         874
#1   10858        3270        5102
#2    5760        3592        1601
#3    4605         577        2925
#4    1911         129          52
#5     551        1217         129
#6     804         298         236
#7     331         639         128

# Create a DataFrame
df <- as.data.frame(cell_counts)

# Sum of each row
row_sums <- rowSums(df)

# Normalize each row
normalized_df <- (df / row_sums)*100

print(normalized_df)

# C19_NDN Healthy_NDN Non_C19_NDN
# 0 27.04068   68.946837    4.012487
# 1 56.46386   17.004680   26.531461
# 2 52.58833   32.794668   14.617000
# 3 56.80276    7.117306   36.079931
# 4 91.34799    6.166348    2.485660
# 5 29.04586   64.153927    6.800211
# 6 60.08969   22.272048   17.638266
# 7 30.14572   58.196721   11.657559

#  Create a DataFrame
df <- as.data.frame(cell_counts)

# Sum of each column
col_sums <- colSums(df)

# Normalize each column
normalized_df <- sweep(df, 2, col_sums, FUN = "/") * 100

print(normalized_df)

# C19_NDN Healthy_NDN Non_C19_NDN
# 0 19.179420  60.7033145    7.911650
# 1 35.356561  13.2174616   46.184484
# 2 18.756106  14.5189976   14.492622
# 3 14.995116   2.3322555   26.477777
# 4  6.222729   0.5214228    0.470716
# 5  1.794204   4.9191593    1.167738
# 6  2.618040   1.2045271    2.136327
# 7  1.077825   2.5828618    1.158686

# Identify clusters and groups
clusters <- levels(covid.filtered.seurat.obj$RNA_snn_res.0.3)
samples <- unique(covid.filtered.seurat.obj$sample)

# Initialize a matrix to store the cell counts
cell_counts <- matrix(0, nrow = length(clusters), ncol = length(samples),
                      dimnames = list(clusters, samples))

# Count the cells in each cluster per group
for (cluster in clusters) {
  for (sample in samples) {
    cell_counts[cluster, sample] <- sum(covid.filtered.seurat.obj$RNA_snn_res.0.3 == cluster &
                                         covid.filtered.seurat.obj$sample == sample)
  }
}

print(cell_counts)

# C19_P1_NDN C19_P2_NDN C19_P3_NDN C19_P4_NDN C19_P5_NDN C19_P6_NDN C19_P7_NDN
# 0        230       2181        304       1777       1398          0          0
# 1        243       3826       1273       3624       1846         26         20
# 2         14       2998       1790        718        197          0         43
# 3        145       2996        428        927        102          0          7
# 4       1017        555        252         69         14          3          1
# 5         37        101         18        194        196          3          2
# 6         29        331        117        164        155          3          5
# 7         26        143         24         83         27         10         18
# Healthy_d1_NDN Healthy_d2_NDN Healthy_d3_NDN Healthy_d4_NDN Healthy_d5_NDN
# 0           7849           7032             40             56             41
# 1           1798           1469              2              0              1
# 2            543           3035              2              4              8
# 3            393            184              0              0              0
# 4             75             53              0              1              0
# 5            549            652              5              8              3
# 6             82            213              2              1              0
# 7            442            180              4              7              6
# Non_C19_P1_NDN Non_C19_P2_NDN
# 0            861             13
# 1           5043             59
# 2           1598              3
# 3           2921              4
# 4             52              0
# 5            120              9
# 6            229              7
# 7            126              2



#  Create a DataFrame
df <- as.data.frame(cell_counts)

# Sum of each column
col_sums <- colSums(df)

# Normalize each column
normalized_df <- sweep(df, 2, col_sums, FUN = "/") * 100

print(normalized_df)

# C19_P1_NDN C19_P2_NDN C19_P3_NDN C19_P4_NDN C19_P5_NDN C19_P6_NDN C19_P7_NDN
# 0 13.2107984 16.6095499  7.2277699 23.5177343 35.5273189   0.000000   0.000000
# 1 13.9574957 29.1371563 30.2662863 47.9618846 46.9123253  57.777778  20.833333
# 2  0.8041356 22.8314675 42.5582501  9.5023822  5.0063532   0.000000  44.791667
# 3  8.3285468 22.8162364 10.1759391 12.2683960  2.5921220   0.000000   7.291667
# 4 58.4147042  4.2266393  5.9914408  0.9131816  0.3557814   6.666667   1.041667
# 5  2.1252154  0.7691722  0.4279601  2.5674960  4.9809403   6.666667   2.083333
# 6  1.6657094  2.5207524  2.7817404  2.1704606  3.9390089   6.666667   5.208333
# 7  1.4933946  1.0890260  0.5706134  1.0984648  0.6861499  22.222222  18.750000
# Healthy_d1_NDN Healthy_d2_NDN Healthy_d3_NDN Healthy_d4_NDN Healthy_d5_NDN
# 0     66.9081920      54.860353      72.727273      72.727273      69.491525
# 1     15.3269116      11.460446       3.636364       0.000000       1.694915
# 2      4.6287614      23.677641       3.636364       5.194805      13.559322
# 3      3.3500980       1.435481       0.000000       0.000000       0.000000
# 4      0.6393317       0.413481       0.000000       1.298701       0.000000
# 5      4.6799079       5.086597       9.090909      10.389610       5.084746
# 6      0.6990026       1.661726       3.636364       1.298701       0.000000
# 7      3.7677947       1.404275       7.272727       9.090909      10.169492
# Non_C19_P1_NDN Non_C19_P2_NDN
# 0      7.8630137      13.402062
# 1     46.0547945      60.824742
# 2     14.5936073       3.092784
# 3     26.6757991       4.123711
# 4      0.4748858       0.000000
# 5      1.0958904       9.278351
# 6      2.0913242       7.216495
# 7      1.1506849       2.061856



covid.filtered.seurat.obj


saveRDS(covid.filtered.seurat.obj, "./NDN_Analysis/my_rds/covid.filtered.seurat.HarmonyIntegration.obj.EmptyDrops.BTfilter.rds")

