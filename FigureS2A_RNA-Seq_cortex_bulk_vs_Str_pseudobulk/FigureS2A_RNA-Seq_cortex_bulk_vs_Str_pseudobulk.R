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

#load files defining the experimental groups and DESeq comparisons to compute
group_tab0 = read_csv(paste0(in_dir,"/group_tab.csv"))
comparisons_tab = read_csv(paste0(in_dir,"/comparisons.csv"))

#define samples and groups
samples = group_tab0$sample
groups = unique(group_tab0$group)

#load count matrices 
input_counts_1 = read_csv(paste0(in_dir,"/matrix.csv"))
input_counts_2 = read_csv(paste0(in_dir,"/01_pseudo_bulk_counts.csv"))


########


####################################################################
####################################################################
# [1] DESeq2 analysis
####################################################################
####################################################################


####################################################################
# get gene symbols for ENSMUS IDs in count matrix (for input_counts_1 (cortex)) 
# then merge with striatum pseudobulk matrix
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

#extract common genes and merge with input_counts_2
t1 = input_counts_1_named
t2 = input_counts_2
common_genes = intersect(t1$gene, t2$gene)
t1 = t1[t1$gene %in% common_genes,]
t1 = t1[order(t1$gene),] 
t2 = t2[t2$gene %in% common_genes,]
t2 = t2[order(t2$gene),] 

input_counts_comb = cbind(t1, t2[,-1])
############

###############################################
# Run DESeq2
###############################################

#generate coldata table for DEseq

coldata1 = as.data.frame(cbind(group = group_tab0$group), row.names = group_tab0$sample)

#get count matrix from count matrix, re-order according to sample order in group_tab0

t1 = as.matrix(input_counts_comb[,-c(1:2)])
row.names(t1) = input_counts_comb$gene
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

expr_tab = as_tibble(cbind(input_counts_comb[,c(1,2)], raw.counts,pseudo.counts, log_matrix))

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

t1 = deseq_results_expressed
comp = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,sep=".")

fc = t1[,paste0(comp,".log2FoldChange")]
p = t1[,paste0(comp,".padj")]

cb = as.matrix(fc)
cb[] = 0
colnames(cb) = paste0("reg_", comp)
cb[(fc <= -1 & p <= 0.05)] = -1
cb[(fc >= 1 & p <= 0.05)] = 1

### extract regulated gene sets

t2 = cbind(deseq_results_expressed, cb)


l1 = list()

l1[["mat_Cor"]] = t2$gene[t2$reg_Cor_P4.vs.Cor_2m == -1]
l1[["immat_Cor"]] = t2$gene[t2$reg_Cor_P4.vs.Cor_2m == 1]

l1[["mat_Str"]] = t2$gene[t2$reg_Str_P3.vs.Str_3m == -1]
l1[["immat_Str"]] = t2$gene[t2$reg_Str_P3.vs.Str_3m == 1]

l1[["mat_common"]] = intersect(l1[["mat_Cor"]], l1[["mat_Str"]])
l1[["immat_common"]] = intersect(l1[["immat_Cor"]], l1[["immat_Str"]])

comp_gene_list = l1


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

write_csv(comp_gene_tab, path = paste0(out_dir,"/02_reg_genes_cor_vs_str_pseudobulk.csv"))

genes_stat = tibble(comparison = names(comp_gene_list), N_genes = lengths(comp_gene_list))
write_csv(genes_stat, path = paste0(out_dir,"/02_reg_genes_cor_vs_str_pseudobulk_stat.csv"))
#############



####################################################################
####################################################################
# [3] downstream analyses/presentations
####################################################################
####################################################################


####################################################################
# extract normalised counts for cell type markers
####################################################################

GOI_tab = read_csv(file = paste0(in_dir,"/Gene of interest table 21-01-20.csv"))

t1 = deseq_results_expressed
v1 = paste0(group_tab0$sample,".pseudo")
t2 = t1[,c("gene", v1)]
t3 = t2[match(na.omit(GOI_tab$cell_type_markers_ext), t2$gene),]

write_csv(t3, file = paste0(out_dir,"/03_pseudocounts_cor_vs_str_cell type markers.csv"))

#########

