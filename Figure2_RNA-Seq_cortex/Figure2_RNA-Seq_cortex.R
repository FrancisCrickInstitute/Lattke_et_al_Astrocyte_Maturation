####################################################
# global settings
####################################################

#define working directories
in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Load packages
library("tidyverse")
library(DESeq2)
library("biomaRt")
library("pheatmap")
library("RColorBrewer")
library(colorRamps)
library(ggrepel)

#load files defining the experimental groups and DESeq comparisons to compute
group_tab0 = read_csv(paste0(in_dir,"/group_tab.csv"))
comparisons_tab = read_csv(paste0(in_dir,"/comparisons.csv"))

#define samples and groups
samples = group_tab0$sample
groups = unique(group_tab0$group)

#load count matrices 
input_counts_1 = read_csv(paste0(in_dir,"/matrix.csv"))
########


####################################################################
####################################################################
# [1] DESeq2 analysis
####################################################################
####################################################################


####################################################################
# get gene symbols for ENSMUS IDs in count matrix
####################################################################

# create expression table base (table with ENSMUS id and gene name)
#   convert ENSMUS identifiers to gene symbols (get translation from biomaRt)

#extract gene identifiers
v1 = input_counts_1$X1

#get corresponding gene symbols from ensembl/bioMart
mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
t2 = mart@attributes
t2 = getBM(attributes=c('ensembl_gene_id','mgi_symbol'), mart = mart)
t3 = t2[match(v1, t2$ensembl_gene_id),]
#replace missing gene symbols by ENSMUS id, duplicated gene symbols by "[symbol]_[ID]"
t3[is.na(t3$mgi_symbol),] = v1[is.na(t3$mgi_symbol)]
t3[t3$mgi_symbol == "",] =  v1[t3$mgi_symbol == ""]
t3[duplicated(t3$mgi_symbol),] = paste(t3$mgi_symbol[duplicated(t3$mgi_symbol)],
                                       t3$ensembl_gene_id[duplicated(t3$mgi_symbol)], 
                                       sep="_")
names(t3)[2] = "gene"

gene_names_tab = t3

# add gene symbols to count matrix
input_counts_1_named = cbind(gene_names_tab, round(input_counts_1[,-1]))

############

###############################################
# Run DESeq2
###############################################

#generate coldata table for DEseq

coldata1 = as.data.frame(cbind(group = group_tab0$group), row.names = group_tab0$sample)

#get count matrix from count matrix, re-order according to sample order in group_tab0

t1 = as.matrix(input_counts_1_named[,-c(1:2)])
row.names(t1) = input_counts_1_named$gene
cts1 = t1[,match(group_tab0$sample, colnames(t1))]

#run DESeq

dds = DESeqDataSetFromMatrix(countData = cts1,
                             colData = coldata1,
                             design = ~ group)
dds = DESeq(dds)

save(dds, file = paste0(out_dir,"/01_Deseq_dataset.rda"))
#######

###############################################
# Extract and filter data
###############################################

#extract raw counts and normalized counts (.pseudo)

raw.counts = counts(dds,normalized=FALSE)
colnames(raw.counts) = paste(colnames(raw.counts),'raw',sep='.')

pseudo.counts = counts(dds,normalized=TRUE)
colnames(pseudo.counts) = paste(colnames(pseudo.counts),'pseudo',sep='.')

# get log2-transformed expression data 

ntd = normTransform(dds)
log_matrix =  as_tibble(assay(ntd))
names(log_matrix) = paste(colnames(log_matrix),'log',sep='.')

# create combined results table

expr_tab = as_tibble(cbind(input_counts_1_named[,c(1,2)], raw.counts,pseudo.counts, log_matrix))

#extract statistics for specified comparisons (in comparisons.csv), add to expr_tab

deseq_results_tab = expr_tab
i=4
for (i in 1:nrow(comparisons_tab)){
  t2 = results(dds, contrast = c("group", comparisons_tab$group1[i], 
                                 comparisons_tab$group2[i]))
  t3 = as_tibble(t2@listData)
  names(t3) = paste(comparisons_tab$group1[i],"vs",comparisons_tab$group2[i],names(t3),sep=".")
  deseq_results_tab = cbind(deseq_results_tab, t3)
}

deseq_results_tab = deseq_results_tab[order(deseq_results_tab$gene),]

write_csv(deseq_results_tab, file = paste0(out_dir,"/01_Deseq_results.csv"))

# remove not expressed genes (max pseudo-expr <10)

t1 = deseq_results_tab
gene_expr_bool = apply(t1[,paste(samples,'pseudo',sep='.')] >=10, 1, any)
deseq_results_expressed = t1[gene_expr_bool,]

write_csv(deseq_results_expressed, path = paste0(out_dir,"/01_Deseq_results_expressed.csv"))

#keep only genes differentially expressed in any pairwise comparison (padj < 0.05)

t1 = deseq_results_expressed
comparisons_padj_names = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,"padj",sep=".")
comparisons_log2FC_names = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,"log2FoldChange",sep=".")
comp_sign_diff = as_tibble((t1[,comparisons_padj_names] <0.05)&(abs(t1[,comparisons_log2FC_names]) > 1))
t2 = t1[apply(comp_sign_diff,1,any)&!apply(comp_sign_diff,1,anyNA),]
t3 = t1[apply(comp_sign_diff,1,anyNA),]

deseq_results_diff = t2

write_csv(deseq_results_diff, path = paste0(out_dir,"/01_Deseq_results_diff_expr.csv"))


#save analysis stats
t1 = tibble(genes = c("all analysed", "expressed", "differentially expressed"),
            N = c(nrow(deseq_results_tab), nrow(deseq_results_expressed), nrow(deseq_results_diff)))
write_csv(t1, path = paste0(out_dir,"/01_Deseq_stats.csv"))
##########



####################################################################
####################################################################
# [2] identify genes regulated by maturation in cortex
####################################################################
####################################################################

####################################################################
# extract regulated genes from DESeq dataset
####################################################################

#create discrete matrix of comparisons up(-1)/down(+1)/not reg(0) for gene list filtering

t1 = deseq_results_diff
comp = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,sep=".")

fc = t1[,paste0(comp,".log2FoldChange")]
p = t1[,paste0(comp,".padj")]

cb = as.matrix(fc)
cb[] = 0
colnames(cb) = paste0("reg_", comp)
cb[(fc <= -1 & p <= 0.05)] = -1
cb[(fc >= 1 & p <= 0.05)] = 1

### extract regulated gene sets

t2 = cbind(deseq_results_diff, cb)

comp_gene_list = NULL

comp_gene_list[["mature_Cor"]] = t2$gene[t2$reg_Astr_P4.vs.Astr_2m < 0]
comp_gene_list[["immature_Cor"]] = t2$gene[t2$reg_Astr_P4.vs.Astr_2m > 0]
#####

####################################################################
# create tables with regulated genes and number of genes in each comparison
####################################################################

l1 = comp_gene_list

comp_gene_tab = as_tibble(matrix(data ="", nrow = max(lengths(l1)), ncol = length(l1), 
                            dimnames = list(NULL, names(l1))))

for (i in names(l1)){
  v1 = l1[[i]]
  comp_gene_tab[c(1:length(v1)),i] = v1
}

write_csv(comp_gene_tab, path = paste0(out_dir,"/02_reg_genes_cortex.csv"))

genes_stat = tibble(comparison = names(comp_gene_list), N_genes = lengths(comp_gene_list))
write_csv(genes_stat, path = paste0(out_dir,"/02_reg_genes_cortex_stat.csv"))
#############



####################################################################
####################################################################
# [3] downstream analyses/presentations
####################################################################
####################################################################

####################################################################
# create volcano plot for comparison with striatal maturation regulated genes
####################################################################

#load regulated genes from striatal sc-RNA-Seq data
comps_str = read_csv(paste0(in_dir,"/06_genes_by_module_cons.csv"))

str_immat_genes = comps_str$gene[comps_str$module_cons %in% c("immat_1", "immat_2", "immat_3", "immat_4")]
str_mat_genes = comps_str$gene[comps_str$module_cons %in% c("mat_1", "mat_2")]

#load genes of interest for labelling on plot

GOI_tab = read_csv(paste0(in_dir,"/Gene of interest table 21-01-20.csv"))

#extract example genes to label
labs_up = na.omit(GOI_tab$mature_sel_volcano)
labs_down = na.omit(GOI_tab$immature_sel_volcano)


#add category(reg in striatum) and labels to deseq data 

t1 = deseq_results_expressed[order(deseq_results_expressed$gene),]

t1$category = "not_reg"
t1$category[t1$gene %in% unlist(comp_gene_list)] = "all_reg_Cor"
t1$category[t1$gene %in% str_mat_genes] = "mature_genes_Str"
t1$category[t1$gene %in% str_immat_genes] = "immature_genes_Str"

v1 = c(labs_down, labs_up)
t1$labels = v1[match(t1$gene, v1)]
t1$labels[is.na(t1$labels)] = ""

deseq_plot_data = t1

#define colors, size and labels for example genes

categories = c("not_reg",  "all_reg_Cor", "mature_genes_Str", "immature_genes_Str")
cat_color = c("grey80", "grey40", "red", "blue")
cat_size = c( 0.2, 0.2, 1, 1)

### volcano plot 

p1 = ggplot(deseq_plot_data, aes(x = -Astr_P4.vs.Astr_2m.log2FoldChange, 
                                 y = -log10(Astr_P4.vs.Astr_2m.padj) ))+
  geom_point(aes(size = category, color = category),alpha = 0.5) + 
  geom_text_repel(aes(label = labels, color = category), fontface = "bold.italic")+
  scale_color_manual(limits = categories, values = cat_color) + 
  scale_size_manual(limits = categories, values = cat_size) + 
  scale_x_continuous(limits = c(-12, 12)) +
  geom_hline(yintercept = -log10(0.05), color = "red", size = 0.3, linetype =2) + 
  geom_vline(xintercept = 1, color = "red", size = 0.3, linetype =2) +
  geom_vline(xintercept = 0, color = "grey60", size = 0.3, linetype =1) +
  geom_vline(xintercept = -1, color = "red", size = 0.3, linetype =2)+ 
  labs(x = "-log2FC", y = "-log10(padj)")+
  theme_bw() +
  theme(text = element_text(size = 14, face =  "bold"),
        axis.text = element_text(size = 14, face = "bold"))


pdf(file = paste0(out_dir,"/03_volcano plot with striatal and selected genes.pdf"), 
    width = 9, height = 5 )
{
  plot(p1)
}

dev.off()

###########

####################################################################
# create heatmap for selected common maturation genes
####################################################################

#load genes of interest and deseq dataset

GOI_tab = read_csv(paste0(in_dir,"/Gene of interest table 21-01-20.csv"))
deseq_results_expressed = read_csv(file = paste0(out_dir,"/01_Deseq_results_expressed.csv"))

# extract log2 expression matrix from deseq table

t1 = deseq_results_expressed
expr_cols = paste0(group_tab0$sample,".log")
expr_matr = as.matrix(t1[,expr_cols]-rowMeans(t1[,expr_cols]))
rownames(expr_matr) = t1$gene

#plot function (to create and save in pdf)

plot_expr_heatmap = function(pl_matr, pl_genes = rownames(pl_matr), 
                             file = "expression_heatmap.pdf",
                             scale_range = NULL, 
                             fontsize = 14, cellwidth = 12, cellheight = 14){
  
  pl_genes = pl_genes[pl_genes %in% rownames(pl_matr)]
  pl_matr = pl_matr[pl_genes,]
  
  if (is.null(scale_range)){
    v1 = max(abs(range(pl_matr, na.rm = TRUE, finite = TRUE)))
    scale_range = c(-v1, +v1)
    }
  
  
  pdf(file = file, width = (3+ncol(pl_matr)*cellwidth/25), 
      height = (3 + nrow(pl_matr))*cellheight/15 )
  {
    pheatmap(pl_matr, cluster_rows=FALSE, show_rownames=TRUE,
             cluster_cols=FALSE, 
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(scale_range[1], scale_range[2], length.out = 251),
             border_color = NA, fontsize = fontsize, fontface = 4,
             cellwidth = cellwidth, cellheight = cellheight,
             main =  paste0("Log2 expression vs mean (", length(pl_genes)," genes not clustered)")
             )
    pheatmap(pl_matr, cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=FALSE, 
             color = colorRampPalette(c("blue", "white", "red"))(250),
             breaks = seq(scale_range[1], scale_range[2], length.out = 251),
             border_color = NA, fontsize = fontsize, fontface = 4,
             cellwidth = cellwidth, cellheight = cellheight,
             main =  paste0("Log2 expression vs mean (", length(pl_genes)," genes clustered)")
             )   
             
  }
  dev.off()
}

#  plot gene sets

plot_expr_heatmap(expr_matr, pl_genes = na.omit(GOI_tab$mature_immature_sel_IF), 
                             file = paste0(out_dir,"/03_heatmap_mature_immature_sel_IF.pdf"))

plot_expr_heatmap(expr_matr, pl_genes = na.omit(GOI_tab$immature_ETS_TFs), scale_range = c(-2,2), 
                  file = paste0(out_dir,"/03_heatmap_mature_immature_sel_ETS_TFs.pdf"))

plot_expr_heatmap(expr_matr, pl_genes = na.omit(GOI_tab$mature_ROR_HOX_TFs), scale_range = c(-2,2), 
                  file = paste0(out_dir,"/03_heatmap_mature_mature_sel_ROR_HOX_TFs.pdf"))


#########

