############################################################################
# global settings
############################################################################

in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

library("tidyverse")
library(DESeq2)
library("biomaRt")
library("pheatmap")

#load count matrix and group table
group_tab = read_csv(paste0(in_dir,"/group_tab.csv"))
comparisons_tab = read_csv(paste0(in_dir,"/comparisons.csv"))

#define samples and groups
samples = group_tab$sample
groups = unique(group_tab$group)
###########

############################################################################
############################################################################
# [1] DESeq2 analysis
############################################################################
############################################################################

############################################################################
# load and format data (count matrices from pipeline)
############################################################################
#load count matrices from RNA-Seq mapping pipeline (round counts to lower integer)
#replace X6 NSC/B14d data with resequenced data (old data were prepared with different kit)

t1 = read_csv(paste0(in_dir,"/matrix_in vitro in vivo_old X6 data.csv"))
t2 = t1
t2[,-1] = floor(t1[,-1])
mat1 = t2

t1 = read_csv(paste0(in_dir,"/matrix_FGF_3D.csv"))
t2 = t1
t2[,-1] = floor(t1[,-1])
mat2 = t2

count_tab_merged = cbind(mat1[,-1], mat2[,-1])
rownames(count_tab_merged) = as_vector(mat1[,1])
##############

############################################################################
# DESeq2 analysis
############################################################################

#generate coldata table for DEseq

coldata1 = as.data.frame(cbind(group = group_tab$group), row.names = group_tab$sample)

#get count matrix from original DESeq dataset, re-order according to sample order in group_tab0
t1 = count_tab_merged
cts1 = t1[,match(group_tab$sample, colnames(t1))]

#run DESeq
dds = DESeqDataSetFromMatrix(countData = cts1,
                             colData = coldata1,
                             design = ~ group)
dds = DESeq(dds)

save(dds, file = paste0(out_dir,"/01_Deseq_dataset.rda"))
#################


############################################################################
# convert ENSMUS IDs to gene symbol, extract normalised counts and DESeq2 statistics
############################################################################

# create expression table base (table with ENSMUS id and gene name)
#   convert ENSMUS identifiers to gene symbols (get translation from biomaRt)

    #extract gene identifiers
v1 = names(dds@rowRanges)

    #for mouse: Mart dataset "mmusculus_gene_ensembl", attribute "mgi_symbol" instead of 'hgnc_symbol'
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


#extract raw counts and normalized counts (.pseudo)

raw.counts = counts(dds,normalized=FALSE)
colnames(raw.counts) = paste(colnames(raw.counts),'raw',sep='.')

pseudo.counts = counts(dds,normalized=TRUE)
colnames(pseudo.counts) = paste(colnames(pseudo.counts),'pseudo',sep='.')


# get log2-transformed expression data 

ntd = normTransform(dds)
log_matrix =  as_tibble(assay(ntd))
names(log_matrix) = paste(colnames(log_matrix),'log',sep='.')

expr_tab = as_tibble(cbind(gene_names_tab, raw.counts,pseudo.counts, log_matrix))


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

write_csv(deseq_results_tab, file = paste0(out_dir,"/01_deseq_results.csv"))
##################

############################################################################
# filter deseq_results to get expressed and differentially expressed genes
############################################################################

# remove not expressed genes (max pseudo-expr <10)

t1 = deseq_results_tab
gene_expr_bool = apply(t1[,paste(samples,'pseudo',sep='.')] >=10, 1, any)
deseq_results_expressed = t1[gene_expr_bool,]

write_csv(deseq_results_expressed, path = paste0(out_dir,"/01_deseq_results_expressed.csv"))


#keep only genes differentially expressed in any pairwise comparison (padj < 0.05, abs(log2FC > 1))

t1 = deseq_results_expressed
comparisons_padj_names = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,"padj",sep=".")
comparisons_log2FC_names = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,"log2FoldChange",sep=".")
comp_sign_diff = as_tibble((t1[,comparisons_padj_names] <0.05)&(abs(t1[,comparisons_log2FC_names]) > 1))
t2 = t1[apply(comp_sign_diff,1,any)&!apply(comp_sign_diff,1,anyNA),]
t3 = t1[apply(comp_sign_diff,1,anyNA),]

deseq_results_diff = t2

write_csv(deseq_results_diff, path = paste0(out_dir,"/01_deseq_results_diff_expr.csv"))

#save gene number stats
t1 = tibble(genes = c("all analysed", "expressed", "differentially expressed"),
            N = c(nrow(deseq_results_tab), nrow(deseq_results_expressed), nrow(deseq_results_diff)))
write_csv(t1, path = paste0(out_dir,"/01_deseq_stats.csv"))

#################


############################################################################
############################################################################
# [2] Classify gene modules according to expression in vivo vs in vitro
############################################################################
############################################################################

############################################################################
# manually define gene modules
############################################################################

# load differential RNA-Seq data (csv)
deseq_diff0 = read_csv(paste0(out_dir,"/01_deseq_results_diff_expr.csv"))

#create binary matrix of differential expression for comparisons up(-1)/down(+1)/not reg(0) for gene list filtering
comp = paste(comparisons_tab$group1,"vs",comparisons_tab$group2,sep=".")

t1 = deseq_diff0
fc = t1[,paste0(comp,".log2FoldChange")]
p = t1[,paste0(comp,".padj")]

cb = fc
names(cb) = paste0("reg_", comp)
cb[] = 0
cb[(fc <= -1 & p <= 0.05)] = -1
cb[(fc >= 1 & p <= 0.05)] = 1

cb = as_tibble(cbind(gene = t1$gene, cb, stringsAsFactors = FALSE))


#### define gene modules, add to gene_list

# all genes regulated by TFs
l1 = list()

l1[["FGF_2D_up"]] = cb$gene[cb$reg_bas_2D.vs.FGF_2D == -1]
l1[["FGF_2D_down"]] = cb$gene[cb$reg_bas_2D.vs.FGF_2D == 1]
l1[["bas_3D_up"]] = cb$gene[cb$reg_bas_2D.vs.bas_3D == -1]
l1[["bas_3D_down"]] = cb$gene[cb$reg_bas_2D.vs.bas_3D == 1]
l1[["FGF_3D_up"]] = cb$gene[cb$reg_bas_2D.vs.FGF_3D == -1]
l1[["FGF_3D_down"]] = cb$gene[cb$reg_bas_2D.vs.FGF_3D == 1]

reg_sets = names(l1) #required to extract gene numbers for pairwise proportion test

# mature genes low in bas_2D controls regulated by TFs

t1 = read_csv(paste0(in_dir,"/07_Cor vs Str immat mat genes_genes.csv"))

v1 = intersect(na.omit(t1$mat_common), cb$gene[cb$reg_bas_2D.vs.Astr_2m == -1])
l2 = lapply(l1, intersect, v1)
names(l2) = paste0("mature_low_in_bas_2D_",names(l1))
l2[["mature_low_in_bas_2D_reg_by_any_signal"]] = unique(unlist(l2))
l2[["mature_low_in_bas_2D"]] = v1

gene_list = c(l1,l2)

lengths(gene_list)
###########

############################################################################
# create gene number stats, csv table with genes in modules, 
# stats of overlap TF regulated genes in mature_low/immature_high genes vs all genes
############################################################################

l1 = gene_list

# gene number stats
t1 = tibble(module = names(l1), N_genes = lengths(l1))
write_csv(t1, paste0(out_dir,"/02_genes_in_modules_stats.csv"))

# csv table with genes in modules
t2 = as_tibble(matrix(data = "", nrow = max(lengths(l1)), ncol = length(l1), 
                   dimnames = list(NULL, names(l1)) ))
for (i in names(l1)){
  if (length(l1[[i]]) > 0){t2[c(1:length(l1[[i]])),i] = l1[[i]]}
}
write_csv(t2, paste0(out_dir,"/02_genes_in_modules.csv"))

# stats of overlap TF regulated genes in mature_low/immature_high genes vs all genes
deseq_results_expressed = read_csv(paste0(out_dir,"/01_deseq_results_expressed.csv"))

N_mat_low = c(lengths(l1[c("mature_low_in_bas_2D" ,paste0("mature_low_in_bas_2D_",reg_sets))]) )
N_expr = c(nrow(deseq_results_expressed), lengths(l1[reg_sets]))

sink(paste0(out_dir,"/02_genes reg in_mat_low vs all expressed_pairwise comp stat.txt"))
{pairwise.prop.test(N_mat_low, N_expr, p.adjust.method = "BH")}
sink()

#############


############################################################################
# plot heatmaps of selected modules
############################################################################

#get gene of interest table with specific gene sets
GOI_tab = read_csv(paste0(in_dir,"/Gene of interest table 21-02-19.csv"))

#get deseq dataset (expressed genes)
deseq_res = read_csv(paste0(out_dir,"/01_deseq_results_expressed.csv"))

### create centered log2 expression matrix from deseq dataset

expr_cols = paste0(group_tab$sample,".log")

m1 = as.matrix(deseq_res[,expr_cols])
rownames(m1) = deseq_res$gene
colnames(m1) = group_tab$sample
m2 = m1 - apply(m1, 1, mean)

expr_matr_centr = m2

### heatmap plot function

expr_heatmap = function(matr, pl_genes = NULL, scale_range = NULL,
                        cellwidth = 16, cellheight = 16, 
                        pl_file = "heatmap.pdf",
                        pl_title = pl_file){
  
  if (nrow(m1)>=2){
    
    if (!is.null(pl_genes)){matr = matr[pl_genes,]}
    if (is.null(scale_range)){
      v1 = max(range(matr), na.rm = TRUE, finite = TRUE)
      scale_range = c(-v1, +v1)
      }
    
    pdf(file = pl_file, width = (2 + ncol(matr)/20*cellwidth), 
        height = (2 + nrow(matr)/30*cellheight) )
    
    {
      pheatmap(matr, cluster_rows = TRUE, show_rownames=TRUE,
               cluster_cols=FALSE, 
               color = colorRampPalette(c("blue", "white", "red"))(250),
               breaks = seq(scale_range[1], scale_range[2], length.out = 251),
               border_color = NA, fontsize = cellheight, fontface = 4,
               cellwidth = cellwidth, cellheight = cellheight,
               main =  paste0(pl_title, " ( ", nrow(matr), " genes)")
      )
      
      pheatmap(matr, cluster_rows = FALSE, show_rownames=TRUE,
               cluster_cols=FALSE, 
               color = colorRampPalette(c("blue", "white", "red"))(250),
               breaks = seq(scale_range[1], scale_range[2], length.out = 251),
               border_color = NA, fontsize = cellheight, fontface = 4,
               cellwidth = cellwidth, cellheight = cellheight,
               main =  paste0(pl_title, " (", nrow(matr), " genes)")
      )
      
    }
    dev.off()
  }
}


### plot selected heatmaps

expr_heatmap(expr_matr_centr, pl_genes = na.omit(GOI_tab$mature_sign_reg_sel), 
             scale_range = NULL, cellwidth = 14, cellheight = 14, 
             pl_file = paste0(out_dir,"/02_heatmap_mature_low_in_bas_2D_reg_by_any_signal_selected.pdf"))

expr_heatmap(expr_matr_centr, pl_genes = na.omit(GOI_tab$mature_TFs_sel), 
             scale_range = NULL, cellwidth = 14, cellheight = 14, 
             pl_file = paste0(out_dir,"/02_heatmap_mature_low_in_bas_2D_TFs_selected.pdf"))


########

