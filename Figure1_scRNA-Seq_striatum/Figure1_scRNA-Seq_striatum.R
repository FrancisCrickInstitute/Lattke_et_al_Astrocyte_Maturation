##############################################################################
# global settings and data load
##############################################################################

#define input/output directories
in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Open packages necessary for analysis.
library(tidyverse)
library(colorRamps)
library(RColorBrewer)
library(Seurat)
library(loomR)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)



#load group table with sample and group information

gr_tab = read_csv(paste0(in_dir,"/group_tab.csv"))
gr = gr_tab$sample

#load genes-of-interest sets

GOI_tab = read_csv("./00_input/Gene of interest table 21-02-19.csv")
##########




#################################################################################################
#################################################################################################
# [1] data load and basic cell clustering
#
# Load single cell expression data and metadata into R
# generate pseudobulk dataset
# generate and filter combined seurat object (full dataset)
#################################################################################################
#################################################################################################

#################################################################################################
# Load single cell expression data and metadata into R
#################################################################################################

# Data derived from 10x can be directly imported into seurat using the Read10X() command followed by path to folder containing output.   
# load all datasets in list, add sample IDs to cell IDs

files_10X = gr_tab$dataset_folder

matrix_list = vector(mode = "list", length = nrow(gr_tab))

for (i in 1:nrow(gr_tab)){
  t1 = Read10X(data.dir = files_10X[i])
  colnames(t1) = paste0(gr_tab$sample[i], "_",colnames(t1))
  matrix_list[[i]] = t1
}

names(matrix_list) = gr_tab$sample
################################

##########################################################################################################
# Create Pseudobulk count table
##########################################################################################################

t1 = NULL

for (i in gr){
  t2 = matrix_list[[i]]
  t3 = apply(t2, 1, sum)
  t1 = cbind(t1, as.numeric(t3))
}

pseudobulk = as_tibble(cbind(names(t3), t1))
names(pseudobulk) = c("gene", gr)

write_csv(pseudobulk, file = paste0(out_dir,"/01_pseudo_bulk_counts.csv"))
###################

##########################################################################################################
# Create Seurat Objects for individual samples
##########################################################################################################

#create seurat objects and add batch metadata
seurat_list = matrix_list

for (i in gr){
  t1 = CreateSeuratObject(counts = matrix_list[[i]], project = gr[i])
  #add metadata from gr_tab
  for (j in 3:ncol(gr_tab)){
    t1[[names(gr_tab)[j] ]] = as.character(gr_tab[gr_tab$sample == i,j])
  }
  t1$sample = i
  seurat_list[[i]] = t1
}
##########

##########################################################################################################
# merge, save full count matrix, metadata, QC metrics and apply QC Filters 
##########################################################################################################

# set thresholds for nGenes and %mito here. 
mitoHi <- 10
nGeneLo <- 500

#merge and calculate mitochondrial content 

seur = merge(x = seurat_list[[1]], y = seurat_list[-1])

seur[["percentMito"]] <- PercentageFeatureSet(seur, pattern = "^mt-")

#save raw counts and metadata
t1 = as.matrix(seur@assays$RNA@counts)
write.csv(t1, paste0(out_dir,"/01_raw_counts_merged.csv"))
t2 = seur@meta.data
write.csv(t2, paste0(out_dir,"/01_raw_counts_merged_cell_metadata.csv"))

#stats for raw dataset by sample 
t2 = read_csv(paste0(out_dir,"/01_raw_counts_merged_cell_metadata.csv"))
t3 = t2 %>% group_by(sample) %>% summarise(N_cells = n(), median(nCount_RNA), 
                                           median(nFeature_RNA), 
                                           N_cells_low_nFeature_RNA = length(nFeature_RNA[nFeature_RNA <= nGeneLo]),
                                           N_cells_high_percentMito = length(nFeature_RNA[percentMito >= mitoHi]),
                                           N_cells_retained = length(nFeature_RNA[nFeature_RNA > nGeneLo & percentMito < mitoHi])
)
write_csv(t3, paste0(out_dir,"/01_raw_merged_dataset_stats_by_sample.csv"))

#filter low quality cells
seur = subset(seur, subset = nFeature_RNA > nGeneLo & percentMito < mitoHi)

save(seur, file = paste0(out_dir,"/01_seurat_filtered.rda"))


#clean memory
rm(t1, t2, matrix_list, seurat_list)
##############

###################################################################################
# normalise expression, run PCA/UMAP and cluster cells, plot UMAP clustering plots (full dataset)
###################################################################################

seur = SCTransform(seur)

seur  <- RunPCA(seur)
seur  <- RunUMAP(seur , dims = 1:30)

seur  <- FindNeighbors(seur , dims = 1:30)
seur  <- FindClusters(seur, verbose = FALSE, resolution = 0.6)


#plot clusters and contribution by sample

gr_colors = brewer.pal(12, name = "Paired")
clust_colors = matlab.like(length(unique(seur$seurat_clusters)))

### plot clusters and contribution by sample (combined)

p1 = DimPlot(seur, reduction = "umap",label = TRUE, pt.size = 0.05, 
             cols = clust_colors) + NoLegend()
p1
dev.off()

p1 = DimPlot(seur, reduction = "umap",label = TRUE, pt.size = 0.05, 
             cols = clust_colors) + NoLegend()

p2 = DimPlot(seur, reduction = "umap", group.by = "group", cols = gr_colors, 
             pt.size = 0.05)+ NoLegend()
p3 = DimPlot(seur, reduction = "umap", group.by = "sample", cols = gr_colors, 
             pt.size = 0.05)

pdf(file = paste0(out_dir,"/01_UMAP_by_cluster_and_sample.pdf"), width = 4, height = 3)
  plot(p1)
  plot(p2)
  plot(p3)
dev.off()
####################

##########################################################################################################
# plot cell type marker expression (full dataset)
##########################################################################################################

p1 = FeaturePlot(seur, features = GOI_tab$cell_type_markers_ext, pt.size = 0.1, 
                 ncol = 6)

pdf(file = paste0(out_dir,"/01_UMAP_cell markers.pdf"), width = 24, height = 18)
  plot(p1)
dev.off()


p1 = DotPlot(seur, features = GOI_tab$cell_type_markers_ext, scale.by = "size") +RotatedAxis()

pdf(file = paste0(out_dir,"/01_Dotplot_cell markers.pdf"), width = 15, height = 8)
  plot(p1)
dev.off()


### plot expression of selected genes on UMAP plots

p1 = FeaturePlot(seur, features = "Sox9", pt.size = 0.05)
p2 = FeaturePlot(seur, features = "Cldn5", pt.size = 0.05)
p3 = FeaturePlot(seur, features = "Sox10", pt.size = 0.05)

pdf(file = paste0(out_dir,"/01_UMAP_Sox9_Cldn5_Sox10.pdf"), width = 4, height = 3)
  plot(p1)
  plot(p2)
  plot(p3)
dev.off()
##################

##########################################################################################################
# save combined dataset, table for cluster labeling, table for selected cluster comparisons 
# (full dataset)
##########################################################################################################

save(seur, file = paste0(out_dir,"/01_seurat_filtered_norm.rda"))

v1 = unique(seur$seurat_clusters)
v2 = v1[order(v1)]
t1 = tibble(cluster = v2, cluster_name = v2)

if (file.exists(paste0(out_dir,"/01_cluster_assignment.csv"))){
  write_csv(t1, paste0(out_dir,"/01_cluster_assignment_new.csv"))
} else {write_csv(t1, paste0(out_dir,"/01_cluster_assignment.csv"))}
############






##########################################################################################################
##########################################################################################################
# [2] cluster labelling (full dataset)
#
# !!! before continuing, assign cluster_name in 01_cluster_assignment.csv !!!
#
# label clusters, create cluster metrics and plots for presentation 
# create table to select clusters for further analysis (02_clusters for reclustering.csv)
##########################################################################################################
##########################################################################################################



##########################################################################################################
# add cluster names, generate cluster stats, 
# plot heatmap and dotplot cell type markers with labeled clusters
# and generate csv for selecting clusters for detailed analysis
##########################################################################################################

clust_assign_tab = read_csv(paste0(out_dir,"/01_cluster_assignment.csv"))

load(file = paste0(out_dir,"/01_seurat_filtered_norm.rda"))


seur$cluster_name = clust_assign_tab$cluster_name[match(seur$seurat_clusters, clust_assign_tab$cluster)]

t1 = seur@meta.data

t3 = t1 %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
N_sample = t1 %>% group_by(sample) %>% summarize(N_cells = n())
t3$N_sample = N_sample$N_cells[match(t3$sample, N_sample$sample)]
t3$fract_sample = t3$N_cells/t3$N_sample
stat_by_sample = t3

### cluster frequency by sample in columns (e.g. for Prism)

freq_by_sample = as_tibble(matrix(nrow = nrow(clust_assign_tab), 
                                  ncol = length(gr)+1))
names(freq_by_sample) = c("Cluster", gr)
freq_by_sample$Cluster = clust_assign_tab$cluster_name
cells_by_sample = freq_by_sample

for (i in gr){
  t1 = stat_by_sample[stat_by_sample$sample == i,]
  freq_by_sample[[i]] = t1$fract_sample[match(freq_by_sample$Cluster,t1$cluster_name )] 
  cells_by_sample[[i]] = t1$N_cells[match(freq_by_sample$Cluster,t1$cluster_name )] 
}

freq_by_sample[is.na(freq_by_sample)] = 0
cells_by_sample[is.na(cells_by_sample)] = 0

write_csv(freq_by_sample, paste0(out_dir,"/02_cell_fract_by_sample.csv"))
write_csv(cells_by_sample, paste0(out_dir,"/02_cell_number_by_sample.csv"))

# plot heatmap and dotplot cell type markers with labeled clusters

p1 = DotPlot(seur, features = GOI_tab$cell_type_markers_ext, group.by = "cluster_name",
             scale.by = "size") +RotatedAxis() + theme(axis.text.x = element_text(face = "italic"))

pdf(file = paste0(out_dir,"/02_Dotplot_cell markers_clusters_named.pdf"), width = 16, height = 8)
plot(p1)
dev.off()

# save seurat with added cluster names
save(seur, file = paste0(out_dir,"/02_seurat_norm_labelled.rda"))

# create csv to select clusters for reclustering

if (file.exists(paste0(out_dir,"/02_clusters for reclustering.csv"))){
  write_csv(clust_assign_tab, paste0(out_dir,"/02_clusters for reclustering_new.csv"))
} else {write_csv(clust_assign_tab, paste0(out_dir,"/02_clusters for reclustering.csv"))}
##################


##########################################################################################################
##########################################################################################################
# [3] Reclustering of selected clusters
#
# !!! delete clusters in '02_clusters for reclustering.csv' to remove from further analyses !!!
##########################################################################################################
##########################################################################################################

###################################################################################
# filter selected clusters, re-normalise expression, re-run PCA/UMAP and re-cluster cells (selected clusters)
###################################################################################

reclust_tab = read_csv(paste0(out_dir,"/02_clusters for reclustering.csv"))

load(file = paste0(out_dir,"/02_seurat_norm_labelled.rda"))


seur = subset(seur, subset = cluster_name %in% reclust_tab$cluster_name)

seur = SCTransform(seur)

seur  <- RunPCA(seur)
seur  <- RunUMAP(seur , dims = 1:30)

seur  <- FindNeighbors(seur , dims = 1:30)
seur  <- FindClusters(seur, verbose = FALSE, resolution = 0.7)


#plot clusters and contribution by sample

gr_colors = brewer.pal(12, name = "Paired")
clust_colors = matlab.like(length(unique(seur$seurat_clusters)))

p1 = DimPlot(seur, reduction = "umap",label = TRUE, pt.size = 0.05, 
             cols = clust_colors) + NoLegend()
p1
dev.off()

p1 = DimPlot(seur, reduction = "umap",label = TRUE, pt.size = 0.05, 
             cols = clust_colors)

p2 = DimPlot(seur, reduction = "umap", group.by = "group", cols = gr_colors, 
             pt.size = 0.05)
p3 = DimPlot(seur, reduction = "umap", group.by = "sample", cols = gr_colors, 
             pt.size = 0.05)

pdf(file = paste0(out_dir,"/03_UMAP_by_cluster_sample_reclustered.pdf"), width = 4.5, height = 3)
plot(p1)
plot(p2)
plot(p3)
dev.off()


### plot expression of selected genes on UMAP plots

p1 = FeaturePlot(seur, features = "Mki67", pt.size = 0.05)

pdf(file = paste0(out_dir,"/03_UMAP_reclustered_Ki67.pdf"), width = 4, height = 3)
plot(p1)
dev.off()


# save seurat of selected clusters reclustered
save(seur, file = paste0(out_dir,"/03_seurat_selected clusters reclustered.rda"))
########


##########################################################################################################
##########################################################################################################
# [4] slingshot pseudotime analysis for identification of lineages 
#
# identification of lineages with cluster assignment (all clusters from [3])
##########################################################################################################
##########################################################################################################

###################################################
# slingshot for all astrocyte clusters 
###################################################

#define clusters for analysis and define color scheme
clusters_for_pseudotime = c(0:16)
cluster_colors = matlab.like(length(clusters_for_pseudotime))

load(file = paste0(out_dir,"/03_seurat_selected clusters reclustered.rda"))


# subset to analyse only main astrocyte clusters

seur = seur[, seur$seurat_clusters %in% clusters_for_pseudotime]

# re-run normalisation and PCA with subset

seur = SCTransform(seur)
seur  <- RunPCA(seur)

#plot PCA with cluster information to identify beginning and end clusters
pdf(file = paste0(out_dir,"/04_Slingshot_PCA astr clusters for lineage constrains.pdf"), width = 7, height = 5)
{
  DimPlot(seur, reduction = "pca", group.by = "seurat_clusters", 
          cols = cluster_colors,
          pt.size = 0.05)
}

dev.off()

#  convert seurat to SingleCellExperiment (as input for Slingshot), and run slingshot

sce = as.SingleCellExperiment(seur)


sce <- slingshot(sce, reducedDim = 'PCA', clusterLabels = sce$seurat_clusters, 
                 start.clus = 6, approx_points = 100)

save(sce, seur, file = paste0(out_dir,"/04_seurat_slingshot analysis_start_clust_6_all_lineages.rda"))
###############

###################################################
#slingshot lineage plots with ggplot2
###################################################

### create plot tables for ggplot lineage curves

# table with median of PC1 and PC2 by cluster

t2 = Embeddings(seur, reduction = "pca")

cluster_PC_median_tab = NULL

for (i in clusters_for_pseudotime){
  t2 = Embeddings(seur, reduction = "pca")[seur$seurat_clusters == i,c(1,2)]
  t3 = tibble(cluster = as.character(i), PC1 = median(t2[,1]), PC2 = median(t2[,2]) )
  cluster_PC_median_tab = rbind(cluster_PC_median_tab, t3)
}

# table with lineage connections of median of PC1 and PC2 by cluster
lineage_mapping_tab = NULL

l1 = SlingshotDataSet(sce)@lineages

for (i in names(l1)){
  for (j in 1:(length(l1[[i]])-1)){
    cluster_start = l1[[i]][j]
    cluster_end = l1[[i]][j+1]
    t2 = tibble(lineage = i, cluster_start = cluster_start, clust_end = cluster_end,
                PC1_start = cluster_PC_median_tab$PC1[cluster_PC_median_tab$cluster == cluster_start], 
                PC2_start = cluster_PC_median_tab$PC2[cluster_PC_median_tab$cluster == cluster_start],
                PC1_end = cluster_PC_median_tab$PC1[cluster_PC_median_tab$cluster == cluster_end], 
                PC2_end = cluster_PC_median_tab$PC2[cluster_PC_median_tab$cluster == cluster_end]
    )
    lineage_mapping_tab = rbind(lineage_mapping_tab, t2)
  }
}

#plot lineages (ggplot)

p1 = DimPlot(seur, reduction = "pca", group.by = "seurat_clusters")+
  geom_segment(data = lineage_mapping_tab, aes(x = PC1_start, y = PC2_start, 
                                              xend = PC1_end, yend = PC2_end), 
              color = "black", size = 0.5)+
  geom_point(data = cluster_PC_median_tab, aes(x = PC1, y = PC2, fill = cluster), 
             size = 1.5, shape = 21)+
  scale_fill_manual(limits = as.character(clusters_for_pseudotime), values = cluster_colors )+
  scale_color_manual(limits = as.character(clusters_for_pseudotime), values = cluster_colors )+
  theme(text = element_text(face = 2, size = 14), 
        axis.text = element_text(face = 2, size = 14))

pdf(file = paste0(out_dir,"/04_seurat_slingshot analysis_start_clust_6_lineages.pdf"), 
    width = 6, height = 4)
  plot(p1)
dev.off()

#create table with clusters in lineages

l1 = SlingshotDataSet(sce)@lineages

t1 = as_tibble(matrix(data = "", nrow = max(lengths(l1)), ncol = length(l1), 
                      dimnames =  list(NULL ,names(l1) ) ))

for (i in names(l1)){
  t1[1:length(l1[[i]]),i] = l1[[i]] 
}

write_csv(t1, paste0(out_dir, "/04_slingshot lineages_clusters included.csv"))
##################




###################################################################################
###################################################################################
# [5] merge selected clusters from [4] with reference dataset from Zeisel et al., 2018 
# (for further cluster characterisation) 
#
###################################################################################
###################################################################################



###################################################################################
# merge selected clusters with reference dataset from Zeisel et al., 2018 
###################################################################################

### prepare own dataset
load(file = paste0(out_dir,"/04_seurat_slingshot analysis_start_clust_6_all_lineages.rda"))

seur_own_data = seur
t1 = seur_own_data@meta.data
seur_own_data$cluster_name = paste0("Str_", seur_own_data$seurat_clusters)
seur_own_data$dataset = "this_study"


### import reference dataset

loom_import <- connect(filename = paste0(in_dir,"/l6_r3_astroependymal_cells.loom"), 
                       mode = "r+")

seur_ref = as.Seurat(loom_import)

loom_import$close_all()

#keep only astr
t1 = seur_ref@meta.data
t2 = t1[t1$TaxonomyRank4 == "Astrocytes",]
seur_ref = subset(seur_ref, cells = rownames(t2))
seur_ref$dataset = "Zeisel_et_al_2018"


### keep only common genes in both sets
common_features = intersect(rownames(seur_own_data), rownames(seur_ref))

seur_own_data = subset(seur_own_data, features = common_features)
seur_ref = subset(seur_ref, features = common_features)

### integrate datasets

seurat_list = list(seur_own_data = seur_own_data, seur_ref = seur_ref)

rm(seur_own_data, seur_ref)

for (i in names(seurat_list)){
  seurat_list[[i]] = SCTransform(seurat_list[[i]])
}

integ.features = SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list = PrepSCTIntegration(object.list = seurat_list, anchor.features = integ.features)

seur_anchors = FindIntegrationAnchors(seurat_list, normalization.method = "SCT", anchor.features = integ.features)

seur_with_ref = IntegrateData(anchorset = seur_anchors, normalization.method = "SCT")

save(seur_with_ref, file = paste0(out_dir,"/05_seurat_with_Zeisel18_astr.rda"))


#remove individual datasets to free memory
rm(seurat_list)
rm(seur_anchors)
##########

###################################################################################
# run PCA/UMAP and cluster cells, save combined seurat object (comparison with Zeisel et al., 2018)
###################################################################################

DefaultAssay(seur_with_ref) <- "integrated"

seur_with_ref  <- RunPCA(seur_with_ref , verbose = FALSE)
seur_with_ref  <- RunUMAP(seur_with_ref , dims = 1:30, verbose = FALSE)

seur_with_ref  <- FindNeighbors(seur_with_ref , dims = 1:30, verbose = FALSE)
seur_with_ref  <- FindClusters(seur_with_ref, verbose = FALSE, resolution = 0.8)

save(seur_with_ref, file = paste0(out_dir,"/05_seurat_with_Zeisel18_astr.rda"))
##########

##########################################################################################################
# plot cluster contributions (comparison with Zeisel et al., 2018)
# include plot with clusters manually merged according to pseudotime lineages
##########################################################################################################

t1 = seur_with_ref@meta.data

t1$orig_cluster = as.numeric(str_remove(t1$cluster_name, "Str_"))
t1$clust_manual = t1$orig_cluster
t1$clust_manual[t1$orig_cluster %in% c(6,8)] = "Str_6_8"
t1$clust_manual[t1$orig_cluster %in% c(0,2)] = "Str_0_2"
t1$clust_manual[t1$orig_cluster %in% c(7, 14)] = "Str_7_14"
t1$clust_manual[t1$orig_cluster %in% c(1,9,11)] = "Str_1_9_11"
t1$clust_manual[t1$orig_cluster %in% c(3,4,5,10,12,15)] = "Str_3_4_5_10_12_15"
t1$clust_manual[t1$orig_cluster %in% c(13, 16)] = "Str_13_16"

seur_with_ref@meta.data = t1

p1 = DimPlot(seur_with_ref, reduction = "umap", group.by = "dataset", 
             cols = matlab.like(length(unique(seur_with_ref$dataset))),
             pt.size = 0.05)

p2 = DimPlot(seur_with_ref, reduction = "umap", group.by = "Description", 
             cols = matlab.like(length(unique(seur_with_ref$Description))),
             pt.size = 0.05, na.value = "grey80")

p3 = DimPlot(seur_with_ref, reduction = "umap", group.by = "Description", 
             cols = matlab.like(length(unique(seur_with_ref$Description))),
             pt.size = 0.05, na.value = "grey80")+NoLegend()

p4 = DimPlot(seur_with_ref, reduction = "umap", group.by = "orig_cluster", 
             cols = cluster_colors, pt.size = 0.05, na.value = "grey80")

p5 = DimPlot(seur_with_ref, reduction = "umap", group.by = "clust_manual", 
             cols = matlab.like(length(na.omit(unique(seur_with_ref$clust_manual)))),
             pt.size = 0.05, na.value = "grey80")

p6 = DimPlot(seur_with_ref, reduction = "umap", group.by = "clust_manual", 
             cols = matlab.like(length(na.omit(unique(seur_with_ref$clust_manual)))),
             pt.size = 0.05, na.value = "grey80")+NoLegend()


pdf(file = paste0(out_dir,"/05_UMAP_Str_astr_clusters_vs_Zeisel18.pdf"), width = 7, height = 5)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
plot(p5)
plot(p6)
dev.off()


rm(seur_with_ref)
#############





##########################################################################################################
##########################################################################################################
# [6] pseudotime gene expression analysis
##########################################################################################################
##########################################################################################################


###################################################
# select main lineage for gene expression analyses, plot pseudotime for main lineage
###################################################

load(file = paste0(out_dir,"/04_seurat_slingshot analysis_start_clust_6_all_lineages.rda"))

main_lineage = "Lineage6"
main_lineage_clusters = slingLineages(sce)[["Lineage6"]]

clusters_for_pseudotime = c(0:16)
cluster_colors = matlab.like(length(clusters_for_pseudotime))

sce$sling_pseudotime = sce[[paste0("slingPseudotime_6")]]
seur$sling_pseudotime = sce$sling_pseudotime

#plot single cells colored by cluster/pseudotime, + cluster centers and lineages with main lineage highlighted

p1 = DimPlot(seur, reduction = "pca", group.by = "seurat_clusters")+
  geom_segment(data = lineage_mapping_tab[lineage_mapping_tab$lineage == main_lineage,], 
               aes(x = PC1_start, y = PC2_start, xend = PC1_end, yend = PC2_end), 
               color = "red", size = 1.5)+
  geom_segment(data = lineage_mapping_tab, 
               aes(x = PC1_start, y = PC2_start, xend = PC1_end, yend = PC2_end), 
               color = "black", size = 0.5)+
  geom_point(data = cluster_PC_median_tab, aes(x = PC1, y = PC2, fill = cluster), 
             size = 2, shape = 21)+
  geom_point(data = cluster_PC_median_tab[cluster_PC_median_tab$cluster %in% main_lineage_clusters,], aes(x = PC1, y = PC2, fill = cluster), 
             size = 2, shape = 21, color = "red")+
  scale_fill_manual(limits = as.character(clusters_for_pseudotime), values = cluster_colors )+
  scale_color_manual(limits = as.character(clusters_for_pseudotime), values = cluster_colors )+
  theme(legend.position = "bottom")+
  labs(title = "Clusters")+guides(fill = "none")

p2 = FeaturePlot(seur, reduction = "pca", features = "sling_pseudotime")+
  scale_color_gradientn(colors = matlab.like(50))+
  geom_segment(data = lineage_mapping_tab[lineage_mapping_tab$lineage == main_lineage,], 
               aes(x = PC1_start, y = PC2_start, xend = PC1_end, yend = PC2_end), 
               color = "red", size = 1.5)+
  geom_segment(data = lineage_mapping_tab, 
               aes(x = PC1_start, y = PC2_start, xend = PC1_end, yend = PC2_end), 
               color = "black", size = 0.5)+
  geom_point(data = cluster_PC_median_tab, aes(x = PC1, y = PC2, fill = cluster), 
             size = 2, shape = 21)+
  geom_point(data = cluster_PC_median_tab[cluster_PC_median_tab$cluster %in% main_lineage_clusters,], aes(x = PC1, y = PC2, fill = cluster), 
             size = 2, shape = 21, color = "red")+
  scale_fill_manual(limits = as.character(clusters_for_pseudotime), values = cluster_colors )+
  theme(legend.position = "bottom", 
        legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+
  guides(fill = "none")


pdf(file = paste0(out_dir,"/06_slingshot analysis_start_clust_6_lineage_6_with pseudotime.pdf"), 
    width = 4, height = 4.5)
plot(p1)
plot(p2)
dev.off()

###########


###################################################
# Identify temporally expressed genes (fitGAM)
###################################################

# isolate clusters in main lineage for analysis of gene expression (only cells with pseudotime assigned for main lineage)

seur = seur[,!is.na(seur$sling_pseudotime)]
sce = sce[,!is.na(sce$sling_pseudotime)]

# fit negative binomial GAM (can take long time)
sce <- fitGAM(counts(sce), cellWeights = rep(1, ncol(sce)), pseudotime = sce$sling_pseudotime)

save(seur, sce, 
     file = paste0(out_dir,"/06_slingshot analysis_start_clust_6_lineage_6_fitGAM_results.rda"))

# test for dynamic expression
ATres <- associationTest(sce)

association_test_tab = as_tibble(cbind(gene = rownames(ATres), ATres))

write_csv(association_test_tab, file = paste0(out_dir,"/06_temp_genes_GAM_association_test.csv"))


###########

###################################################
# extract gene expression matrix and bins for expression along pseudotime 
# (for plotting and identification of expression modules)
# 
# define empirically number of bins
# create statistics of contribution of cell clusters to specific pseudotime bins
###################################################

load(file = paste0(out_dir,"/06_slingshot analysis_start_clust_6_lineage_6_fitGAM_results.rda"))
association_test_tab = read_csv(file = paste0(out_dir,"/06_temp_genes_GAM_association_test.csv"))

clusters_for_pseudotime = c(0:16)
cluster_colors = matlab.like(length(clusters_for_pseudotime))

n_bins = 20

### extract SCT norm expression matrix and cell metadata from seurat object

t0 = seur@meta.data
t1 = as_tibble(cbind(cell.id = as.character(rownames(t0)), t0))
cell_meta = t1[order(t1$sling_pseudotime),]

pl_cells = as.character(cell_meta$cell.id)

m1 = seur@assays$SCT@data
m2 = m1[,colnames(m1) %in% pl_cells]
expr_mat = m2[,order(match(colnames(m2), pl_cells))]


### bin cells in same size pseudotime bins and create stats for contribution of seurat clusters to each bin
# create mean expression per bin (only consider bins with >= 10 cells), normalise to mean over all bins

m0 = expr_mat

m1 = as.matrix(m0[rownames(m0) %in% association_test_tab$gene,])

clust_expr_mat = matrix(nrow = nrow(m1), ncol = n_bins, dimnames = list(rownames(m1), 1:n_bins))

max_pseudotime = max(cell_meta$sling_pseudotime)
pseudotime_bin_size = max_pseudotime/n_bins

pseudotime_cluster_stat = NULL
seur$pseudotime_bin = NA_integer_

for (i in 1 : n_bins){
  
  bin_cells = cell_meta$cell.id[(cell_meta$sling_pseudotime > (i-1)*pseudotime_bin_size & 
                                   cell_meta$sling_pseudotime <= i*pseudotime_bin_size)]
  # add cell bin information to seur object
  seur$pseudotime_bin[colnames(seur) %in% bin_cells] = i
  # calculate mean expression per gene in bin
  if (length(bin_cells)>10){
    m2 = m1[,colnames(m1) %in% bin_cells]
    clust_expr_mat[,i] = apply(m2, 1, mean, na.rm = TRUE)
  }
  
  #create stats of cells by cluster and bin
  t1 = cell_meta[cell_meta$cell.id %in% bin_cells,]
  t2 = t1 %>% group_by(seurat_clusters) %>% summarize(n())
  t3 = tibble(cluster = clusters_for_pseudotime, 
              bin = i, bin_max_pseudotime = i*pseudotime_bin_size,
              n_cells = t2$`n()`[match(clusters_for_pseudotime, t2$seurat_clusters)])
  t3$fract_of_bin = t3$n_cells/sum(na.omit(t3$n_cells))
  pseudotime_cluster_stat = rbind(pseudotime_cluster_stat, t3)
}

pseudotime_cluster_stat$cluster = as.character(pseudotime_cluster_stat$cluster)
pseudotime_cluster_stat[is.na(pseudotime_cluster_stat$n_cells), 
                        c("n_cells", "fract_of_bin")] = 0

write_csv(pseudotime_cluster_stat, path = paste0(out_dir,"/06_stat_cluster_by_pseudotime_bin.csv"))


### scale clust_expr_mat and remove lowly changing genes

m1 = clust_expr_mat
#mean-center and remove lowly changing genes
m2 = m1 - apply(m1, 1, mean, na.rm = TRUE)
m3 = m2[apply(abs(m2),1, max, na.rm = TRUE)>0.2,]

clust_expr_mat_cent = m3

#quick check specific genes
plot(clust_expr_mat_cent[rownames(clust_expr_mat_cent) == "Fezf2", ])
########

###################################################
# Identify initial gene modules with similar expression pattern  
# 
# clustering of binned expression matrix with pheatmap function
# define empirically number initial gene modules for manual consolidation
# create files with genes and TFs in specific modules
# create statistics of number of genes in each module
# create file for manual cluster consolidation (module_consolidation_table.csv)
###################################################

N_gene_modules = 50

#create plot for clustering (high intensity to improve clustering for lowly changing genes)

m1 = clust_expr_mat_cent

max_range = max(range(is.finite(m1)))
lim = c(-max_range, max_range )

mod_plot = pheatmap(m1, show_rownames=TRUE, cluster_rows = TRUE,
                    cluster_cols = FALSE, show_colnames = TRUE, 
                    clustering_distance_rows = "euclidean",
                    clustering_method = "ward.D2",
                    treeheight_row = 100,
                    cutree_rows = N_gene_modules,
                    color = colorRampPalette(c("blue", "white", "red"))(250),
                    breaks = seq(lim[1]/5, lim[2]/5, length.out = 251),
                    border_color = NA, fontsize = 10, fontface = 2,
                    cellwidth = 2, cellheight = 0.1
)

dev.off()

#extract plot order and cluster information from module plot object (mod_plot) 
t2 = cutree(mod_plot$tree_row, k=N_gene_modules)
t3 = tibble(gene = names(t2), module = as.character(t2))
t3 = t3[mod_plot$tree_row$order,]
mod_assign_tab = t3


### save plot with cluster assignments

m1 = clust_expr_mat_cent[match(mod_assign_tab$gene, rownames(clust_expr_mat_cent)),]
rownames(m1) = paste0(rownames(m1)," (",mod_assign_tab$module,")")

max_range = max(range(is.finite(m1)))
lim = c(-max_range, max_range )
vert_size = 0.5

pdf(file = paste0(out_dir,"/06_heatmap_by_pseudotime_with_gene_modules.pdf"), width = 16, 
    height = nrow(m1)/80+1)
{
  pheatmap(m1, show_rownames=TRUE, cluster_rows = TRUE,
           cluster_cols = FALSE, show_colnames = TRUE, 
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 100,
           cutree_rows = N_gene_modules,
           color = colorRampPalette(c("blue", "white", "red"))(250),
           breaks = seq(lim[1], lim[2], length.out = 251),
           border_color = NA, fontsize = 5, fontface = 2,
           cellwidth = 10, cellheight = vert_size,
           main =  "expression modules (SCT norm expression, low intens)"
  )
  
  pheatmap(m1, show_rownames=TRUE, cluster_rows = TRUE,
           cluster_cols = FALSE, show_colnames = TRUE, 
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           treeheight_row = 100,
           cutree_rows = N_gene_modules,
           color = colorRampPalette(c("blue", "white", "red"))(250),
           breaks = seq(lim[1]/5, lim[2]/5, length.out = 251),
           border_color = NA, fontsize = 5, fontface = 2,
           cellwidth = 10, cellheight = vert_size,
           main =  "expression modules (SCT norm expression, high intens)"
  )
  
  
}
dev.off()


#save cluster assignment of genes and TFs

write_csv(mod_assign_tab, path = paste0(out_dir,"/06_genes_by_module.csv"))

t1 = mod_assign_tab[mod_assign_tab$gene %in% GOI_tab$TFs_Fantom5_mouse,] 
write_csv(t1, path = paste0(out_dir,"/06_TFs_by_module.csv"))


#generate list of genes by module

mod_list = list()

for (i in unique(mod_assign_tab$module)){
  mod_list[[i]] = mod_assign_tab$gene[mod_assign_tab$module == i]
}


#generate table for manual cluster merging/consolidation

mod_cons_tab = tibble(module = unique(mod_assign_tab$module), 
                      module_consol = unique(mod_assign_tab$module))

if (file.exists(paste0(out_dir,"/06_module_consolidation_table.csv"))){
  write_csv(mod_cons_tab, path = paste0(out_dir,"/06_module_consolidation_table_new.csv"))
} else {write_csv(mod_cons_tab, path = paste0(out_dir,"/06_module_consolidation_table.csv"))}


# save seur with processed expression matrices and gene module assignment

save(seur, sce, expr_mat, clust_expr_mat_cent, 
     mod_list, mod_assign_tab, association_test_tab,
     file = paste0(out_dir,"/06_seurat_slingshot_analysis_with_binned_gene_expression_modules.rda"))

#########

###################################################
# add manual names in module_consolidation_table.csv for module consolidation
# (merge modules with similar expression) 
###############################################

###################################################
# expression module consolidation
#
# add consolidated module to mod_assign_tab
# generate list and stats of genes and TFs by cons module 
###################################################

mod_cons_tab = read_csv(paste0(out_dir,"/06_module_consolidation_table.csv"))

# add consolidated clusters to clust_assign_tab (and save cluster assignment genes/TFs)
t1 = mod_assign_tab

t1$module_cons = mod_cons_tab$module_consol[match(t1$module, mod_cons_tab$module)]
t1 = t1[order(t1$module_cons),]
mod_assign_tab = t1

write_csv(t1, path = paste0(out_dir,"/06_genes_by_module_cons.csv"))

t2 = t1[t1$gene %in% GOI_tab$TFs_Fantom5_mouse,] 

write_csv(t2, path = paste0(out_dir,"/06_TFs_by_module_cons.csv"))


#generate list of genes by cons cluster

mod_list_cons = NULL

for (i in unique(mod_assign_tab$module_cons)){
  mod_list_cons[[i]] = mod_assign_tab$gene[mod_assign_tab$module_cons == i]
}


### generate stats of genes by consolidated module

t2 = mod_assign_tab
t3 = t2 %>% group_by(module_cons) %>% summarise(N_genes = n())
write_csv(t3, path = paste0(out_dir,"/06_Stats_genes_by_module_cons.csv"))


# save modified data set

save(seur, sce, expr_mat, clust_expr_mat_cent, 
     mod_list, mod_assign_tab, association_test_tab, mod_list_cons, mod_cons_tab,  
     file = paste0(out_dir,"/06_seurat_slingshot_analysis_with_binned_gene_expression_modules.rda"))
########

###################################################
# plot presentation heatmap of consolidated modules and graph of cells per bin by cluster
###################################################

load(file = paste0(out_dir,"/06_seurat_slingshot_analysis_with_binned_gene_expression_modules.rda"))
pseudotime_cluster_stat = read_csv(paste0(out_dir,"/06_stat_cluster_by_pseudotime_bin.csv"))
pseudotime_cluster_stat$cluster = as.factor(pseudotime_cluster_stat$cluster)

#define clusters for analysis and define color scheme
clusters_for_pseudotime = c(0:16)
cluster_colors = matlab.like(length(clusters_for_pseudotime))

# remove bins 1 and 20 (too few cells) from expression matrix and cell by bin stats
# reorder modules for plotting

clust_expr_mat_cent_sel = clust_expr_mat_cent[,c(2:19)]

pseudotime_cluster_stat = pseudotime_cluster_stat[pseudotime_cluster_stat$bin %in% c(2:19),]

mod_order = c(
  "immat_1" ,"immat_2", "immat_3" ,"immat_4" , "mat_1"  , "mat_2"
)


### plot cell cluster contribution per bin

p1 = ggplot(data = pseudotime_cluster_stat) + 
  theme(text = element_text(face = 2, size = 28))

p2 = p1 + geom_line(aes(x = bin,y = n_cells, color = cluster ), size = 1.5) +
  scale_color_manual(limits = clusters_for_pseudotime, 
                     values = cluster_colors)+
  scale_x_discrete(limits = pseudotime_cluster_stat$bin)

p3 = p1 + geom_line(aes(x = bin,y = fract_of_bin, color = cluster ), size = 1.5) +
  scale_color_manual(limits = clusters_for_pseudotime, 
                     values = cluster_colors)+
  scale_x_discrete(limits = pseudotime_cluster_stat$bin)


pdf(file = paste0(out_dir,"/06_stat_cluster_by_pseudotime_bin_2-19.pdf"), width = 6, 
    height = 2)
{
  plot(p2)
  plot(p3)
}
dev.off()



### plot gene expression by pseudotime bin; 
# overview for all regulated genes, TFs and selected genes for validation 

plot_module_matrix = function(pl_matr, mod_assign_tab, pl_genes = NULL,
                              plot_file = "expression_matrix.pdf", plot_modules = FALSE,
                              vert_size = 10, font_size = 10, cellwidth = 5){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_matr)}
  
  #reorder and re-label expr matrix
  t1 = mod_assign_tab
  t2 = t1[t1$gene %in% pl_genes,]
  
  m1 = pl_matr[match(t2$gene, rownames(pl_matr)),]
  
  # define scale limits and vertical plot size
  v1 = max(range(m1, na.rm = TRUE, finite = TRUE))
  lim = c(-v1, v1)
  
  ### assign gaps and labels if plot_modules = TRUE
  
  if (plot_modules == TRUE){
    rownames(m1) = paste0(rownames(m1)," (",t2$module_cons,")")
    
    # define gaps at cluster borders
    t3 = t2 %>% group_by(module_cons) %>% summarise(N_genes = n())
    t3 = t3[order(match(t3$module_cons, unique(t2$module_cons))),]
    gaps_row = t3$N_genes[1]
    for (i in 2:(nrow(t3)-1)){gaps_row = c(gaps_row,last(gaps_row) + t3$N_genes[i])}
    
    # define row annotation (cons clusters and original clusters)
    annotation_row = as.data.frame(t2[,-1]) #does not work with tibble
    rownames(annotation_row) = rownames(m1)
  } else {
    gaps_row = NULL
    m1 = m1[order(match(rownames(m1), pl_genes)),]
  }
  
  
  pdf(file = plot_file, width = 8, 
      height = nrow(m1)/30*vert_size+1)
  {
    pheatmap(m1, show_rownames=TRUE, cluster_rows = FALSE,
             cluster_cols = FALSE, show_colnames = TRUE, 
             gaps_row = gaps_row, annotation_row = NULL,
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(lim[1], lim[2], length.out = 251),
             border_color = NA, fontsize = font_size, fontface = 4,
             cellwidth = cellwidth, cellheight = vert_size,
             main =  "expression by pseudotime (SCT norm expr, low intens)"
    )
    
    pheatmap(m1, show_rownames=TRUE, cluster_rows = FALSE,
             cluster_cols = FALSE, show_colnames = TRUE, 
             gaps_row = gaps_row, annotation_row = NULL,
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(lim[1]/5, lim[2]/5, length.out = 251),
             border_color = NA, fontsize = font_size, fontface = 4,
             cellwidth = cellwidth, cellheight = vert_size,
             main =  "expression by pseudotime (SCT norm expr, high intens)"
    )
    
  }
  dev.off()
  
}  


m1 = clust_expr_mat_cent_sel

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab, 
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_cons_modules.pdf"),
                   plot_modules = TRUE,
                   vert_size = 0.1, font_size = 16)

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab,
                   pl_genes = GOI_tab$TFs_Fantom5_mouse,
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_all_TFs.pdf"),
                   plot_modules = TRUE,
                   vert_size = 10, font_size = 10)

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab,
                   pl_genes = GOI_tab$mature_immature_sel_IF,
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_mat_immat_sel_IF.pdf"),
                   vert_size = 10, font_size = 10)

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab,
                   pl_genes = GOI_tab$immature_ETS_TFs,
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_mat_immat_sel_ETS_TFs.pdf"),
                   vert_size = 10, font_size = 10)

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab,
                   pl_genes = GOI_tab$mature_ROR_HOX_TFs,
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_mat_mat_sel_ROR_HOX_TFs.pdf"),
                   vert_size = 10, font_size = 10)

plot_module_matrix(pl_matr = m1, 
                   mod_assign_tab = mod_assign_tab,
                   pl_genes = GOI_tab$mature_TFs_sel,
                   plot_file = paste0(out_dir,"/06_heatmap_pseudotime_gene_expr_mat_mat_TF_sel_BMP_low.pdf"),
                   vert_size = 10, font_size = 10)


#cleanup
rm(m1,m2,m3, clust_expr_mat, clust_expr_mat_cent, clust_expr_mat_cent_sel,mod_plot, 
   expr_mat, seur, sce)


##########



##########################################################################################################
##########################################################################################################
# [7] downstream analysis of immature/mature genes identified by pseudotime analysis
##########################################################################################################
##########################################################################################################


###################################################
# GO analysis of immature and mature genes (clusterprofiler)
###################################################

load(file = paste0(out_dir,"/06_seurat_slingshot_analysis_with_binned_gene_expression_modules.rda"))

# merge different immat/mat gene sets
l1 = mod_list_cons[c("immat_1", "immat_2", "immat_3", "immat_4")]
l2 = mod_list_cons[c("mat_1"  , "mat_2")]

GO_genes_list = list(immature = unlist(l1),
                     mature = unlist(l2))

#perform GO analysis

for (i in names(GO_genes_list)){
  
    ego = enrichGO(gene         = GO_genes_list[[i]],
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

    #save enrichment results
    
    write_csv(ego@result, path = paste0(out_dir,"/07_GO_",i,".csv"))
    
}

#########

###################################################
# overlap of striatal immature/mature genes with cortical genes from bulk RNA-Seq
###################################################

### select comparisons

comp_list = list()

t1 = read_csv("./00_input/02_reg_genes_cortex.csv")
comp_list[["immat_Cor"]] = na.omit(t1$immature_Cor)
comp_list[["mat_Cor"]] = na.omit(t1$mature_Cor)

l1 = mod_list_cons
names(l1)
l2 =l1[c("immat_1", "immat_2", "immat_3", "immat_4")]
comp_list[["immat_Str"]] = unlist(l2)

l2 = l1[c("mat_1" ,  "mat_2")]
comp_list[["mat_Str"]] = unlist(l2)

comp_list[["mat_common"]] = intersect(comp_list[["mat_Str"]], comp_list[["mat_Cor"]])
comp_list[["immat_common"]] = intersect(comp_list[["immat_Str"]], comp_list[["immat_Cor"]])

comp_list[["mat_common_TFs"]] = intersect(comp_list[["mat_common"]], GOI_tab$TFs_Fantom5_mouse)
comp_list[["immat_common_TFs"]] = intersect(comp_list[["immat_common"]], GOI_tab$TFs_Fantom5_mouse)

v1 = comp_list[["immat_Cor"]]
v2 = v1[!(v1 %in% comp_list[["immat_Str"]])]
comp_list[["immat_Cor_only_TFs"]] = intersect(v2, GOI_tab$TFs_Fantom5_mouse)
v1 = comp_list[["mat_Cor"]]
v2 = v1[!(v1 %in% comp_list[["mat_Str"]])]
comp_list[["mat_Cor_only_TFs"]] = intersect(v2, GOI_tab$TFs_Fantom5_mouse)
comp_list[["mat_Cor_only_HOX_TFs"]] = intersect(v2, GOI_tab$Homeo_TFs_mouse)

v1 = comp_list[["immat_Str"]]
v2 = v1[!(v1 %in% comp_list[["immat_Cor"]])]
comp_list[["immat_Str_only_TFs"]] = intersect(v2, GOI_tab$TFs_Fantom5_mouse)
v1 = comp_list[["mat_Str"]]
v2 = v1[!(v1 %in% comp_list[["mat_Cor"]])]
comp_list[["mat_Str_only_TFs"]] = intersect(v2, GOI_tab$TFs_Fantom5_mouse)
comp_list[["mat_Str_only_HOX_TFs"]] = intersect(v2, GOI_tab$Homeo_TFs_mouse)


### generate stats by comparison

t1 = tibble(comp = names(comp_list), N_genes = lengths(comp_list))
write_csv(t1, paste0(out_dir,"/07_Cor vs Str immat mat genes_stats.csv"))


### save genes by comparison as table

l1 = comp_list
t1 = as_tibble(matrix(data = "",nrow = max(lengths(l1)), ncol = length(l1)))
names(t1) = names(l1)

for (i in names(l1)){
  v1 = l1[[i]]
  if (length(v1)>0){
    t1[1:length(v1),i] = v1
  }
}

write_csv(t1, paste0(out_dir,"/07_Cor vs Str immat mat genes_genes.csv"))


########


###################################################
# analysis of expression of immature/mature genes compared with reference dataset 
# from Zeisel et al., 2018
##################################################

# load data and define genes of interest
GOI_tab_2 = read_csv(paste0(out_dir,"/07_Cor vs Str immat mat genes_genes.csv"))

load(file = paste0(out_dir,"/05_seurat_with_Zeisel18_astr.rda"))

#### re-normalise with SCTransform

seur_with_ref = SCTransform(seur_with_ref)

save(seur_with_ref, file = paste0(out_dir,"/07_seurat_with_Zeisel18_astr.rda"))


# define cell clusters for plotting/analysis

load(file = paste0(out_dir,"/07_seurat_with_Zeisel18_astr.rda"))

seur_with_ref$cluster_name_merged = seur_with_ref$cluster_name
seur_with_ref$cluster_name_merged[is.na(seur_with_ref$cluster_name_merged)] =
  seur_with_ref$Description[is.na(seur_with_ref$cluster_name_merged)]

clusters0 = unique(seur_with_ref$cluster_name_merged)
clusters0
clusters = c(paste0("Str_", c(6,8,2,0,7,14,1,3,4,5,9,10,11,12,15,13,16)), 
             "Telencephalon astrocytes, fibrous" ,             "Telencephalon astrocytes, protoplasmic",        
             "Olfactory astrocytes" ,         
             "Non-telencephalon astrocytes, protoplasmic" ,    "Non-telencephalon astrocytes, fibrous",
             "Dorsal midbrain Myoc-expressing astrocyte-like", "Bergmann glia"
)

### extract SCT norm expression matrix and cell metadata from seurat object 

t0 = seur_with_ref@meta.data
t1 = as_tibble(cbind(cell.id = as.character(rownames(t0)), t0))
t2 = t1[t1$cluster_name_merged %in% clusters,]
cell_meta = t2[order(match(t2$cluster_name_merged, clusters)),]
pl_cells = as.character(cell_meta$cell.id)

m1 = seur_with_ref@assays$SCT@data
m2 = m1[,colnames(m1) %in% pl_cells]
expr_mat = m2[,order(match(colnames(m2), pl_cells))]

# calculate mean expr by cluster

m1 = expr_mat
clust_expr_mat = matrix(nrow = nrow(m1), ncol = length(clusters), 
                        dimnames = list(rownames(m1), clusters))

for (i in clusters){
  clust_cells = cell_meta$cell.id[cell_meta$cluster_name_merged == i]
  m2 = m1[,colnames(m1) %in% clust_cells]
  clust_expr_mat[,i] = apply(m2, 1, mean)
}

#mean-center clust_expr_mat
m1 = clust_expr_mat
m2 = m1 - apply(m1, 1, mean)
clust_expr_mat_cent = m2

### plot heatmaps for immature/mature genes

pl_genes = na.omit(GOI_tab_2$mat_common)
pl_matr = clust_expr_mat_cent
rm(pl_genes, pl_matr, lim)

plot_expr_matrix = function(pl_matr, pl_genes = NULL, 
                            plot_file = "heatmap_expression_matrix.pdf", 
                            pl_title = "", lim = NULL,
                              vert_size = 10, font_size = 10, cellwidth = 10){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_matr)}
  
  #reorder expr matrix
  m1 = pl_matr[match(pl_genes, rownames(pl_matr), nomatch = 0),]
  
  # auto-define scale limits and vertical plot size if not given
  if (is.null(lim)){
    v1 = max(abs(range(m1, na.rm = TRUE, finite = TRUE)))
    lim = c(-v1, v1)
  }
  
  pdf(file = plot_file, width = 8, 
      height = nrow(m1)/30*vert_size+1)
  {
    pheatmap(m1, show_rownames=TRUE, cluster_rows = FALSE,
             cluster_cols = FALSE, show_colnames = TRUE, 
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(lim[1], lim[2], length.out = 251),
             border_color = NA, fontsize = font_size, fontface = 4,
             cellwidth = cellwidth, cellheight = vert_size,
             main =  pl_title
    )
    pheatmap(m1, show_rownames=TRUE, cluster_rows = TRUE,
             cluster_cols = TRUE, show_colnames = TRUE, 
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(lim[1], lim[2], length.out = 251),
             border_color = NA, fontsize = font_size, fontface = 4,
             cellwidth = cellwidth, cellheight = vert_size,
             main =  pl_title
    )
    pheatmap(m1, show_rownames=TRUE, cluster_rows = TRUE,
             cluster_cols = FALSE, show_colnames = TRUE, 
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(lim[1], lim[2], length.out = 251),
             border_color = NA, fontsize = font_size, fontface = 2,
             cellwidth = cellwidth, cellheight = vert_size,
             main =  pl_title
    )
  }
  dev.off()
}  


plot_expr_matrix(pl_matr = clust_expr_mat_cent, pl_genes = na.omit(GOI_tab_2$mat_common), 
                  plot_file = paste0(out_dir,"/07_heatmap_expr_vs_Zeisel18_mat_common.pdf"), 
                  vert_size = 0.5, lim = c(-1,1))

plot_expr_matrix(pl_matr = clust_expr_mat_cent, pl_genes = na.omit(GOI_tab_2$immat_common), 
                 plot_file = paste0(out_dir,"/07_heatmap_expr_vs_Zeisel18_immat_common.pdf"), 
                 vert_size = 0.5, lim = c(-1,1))




########


##########################
# end, save session info
##########################
session_info = sessionInfo()
save(session_info, file = paste0(out_dir,"/sessionInfo.rda"))
writeLines(capture.output(session_info), paste0(out_dir,"/sessionInfo.txt"))


