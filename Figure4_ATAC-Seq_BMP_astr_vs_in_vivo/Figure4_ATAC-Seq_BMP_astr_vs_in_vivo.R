##############################################################################
# global settings
##############################################################################

in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(tidyverse)
library(DESeq2)

###############


##############################################################################
##############################################################################
# [1] peak annotations from ATAC dataset with reference datasets 
##############################################################################
##############################################################################


##############################################################################
# load ATAC peak annotations from pipeline and reference datasets for annotations 
##############################################################################

peak_annot0 = read_tsv(paste0(in_dir,"/merged_peaks.homer.annotatePeaks.txt"))

### extract coordinates and annotations from ATAC peak annotations from pipeline 

t1 = peak_annot0
names(t1)[1] = "Peak_id"
v1 = c("Peak_id", "Chr", "Start", "End", "Strand", 
       "Annotation", "Distance to TSS", "Nearest PromoterID")
t2 =  t1[v1]
peak_annot1 = t2[order(t2$Chr, t2$Start, t2$End),]

# define Chromosome annotations for which data will be analysed (discards regions with unknown location and patches)

Chr = paste0("chr",c(1:19, "X", "Y"))


### load bed files with features for additional annotations
#   ensemble regulatory annotations and H3K4me1 enhancer mark ChIP peaks from encode database (forebrain P0 and adult cortex)
#   peaks from brain DNase-Seq from encode
#   enhancer/promoter pairings from Shen12 and Ron17 for enhancer-gene assignment 
#   format for overlap: 4 columns: c("Chr", "Start", "End", "Feature"), add to feat_list
#       "Feature" value will be added to original peak annotation peak_annot1 
#       (column name from element name in feat_list)

feat_list = list()

#extract location and annotation from ensembl regulatory built
t1 = read_tsv(file = "/Users/lattkem/Desktop/Omics analysis/genomic references/Ensembl references/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20161111.gff", col_names = FALSE)
t2 = t1[,c(1,4,5,9)]
names(t2) = c("Chr", "Start", "End", "Descr") 
t2$Chr = paste0("chr",t2$Chr)
t3 = str_split(t2$Descr, ";", simplify = TRUE)
t4 = as_tibble(t3)
v1 = str_remove(t4$V5, "feature_type=")
t5 = cbind(t2[,c(1:3)], Feature = v1)
feat_list[["Ensemble_feature_type"]] = t5

t1 = read_tsv(file = "/Users/lattkem/Desktop/Omics analysis/genomic references/reference ChIP Encode/ENCFF172LKQ_Fo_P0_H3K4me1.bed", col_names = FALSE)
t2 = t1[,c(1:4)]
t2[,4] = "H3K4m1_Fo_P0"
feat_list[["H3K4m1_Fo_P0"]] = t2

t1 = read_tsv(file = "/Users/lattkem/Desktop/Omics analysis/genomic references/reference ChIP Encode/ENCFF746YEV_H3K4me1_Cor_8w.bed", col_names = FALSE)
t2 = t1[,c(1:4)]
t2[,4] = "H3K4m1_Cor_8w"
feat_list[["H3K4m1_Cor_8w"]] = t2

t1 = read_tsv(file = "/Users/lattkem/Desktop/Omics analysis/genomic references/reference bigwigs Encode brain DNase-Seq/ENCFF591XUM_brain E18.bed", col_names = FALSE)
t2 = t1[,c(1:4)]
t2[,4] = "DHS_brain_E18"
feat_list[["DHS_brain_E18"]] = t2

t1 = read_tsv(file = "/Users/lattkem/Desktop/Omics analysis/genomic references/reference bigwigs Encode brain DNase-Seq/ENCFF865BUI_brain_adult.bed", col_names = FALSE)
t2 = t1[,c(1:4)]
t2[,4] = "DHS_brain_ad"
feat_list[["DHS_brain_ad"]] = t2

t1 = read_tsv(file = paste0(in_dir,"/Shen12_brainE14_mm10.bed"),col_names = FALSE, col_types = "cnncc") 
t2 = t1[,c(1:4)]
t2$X1 = paste0("chr",t2$X1)
feat_list[["Enh_Shen12_brain_E14"]] = t2

t1 = read_tsv(file = paste0(in_dir,"/Shen12_cortex_mm10.bed"),col_names = FALSE, col_types = "cnncc") 
t2 = t1[,c(1:4)]
t2$X1 = paste0("chr",t2$X1)
feat_list[["Enh_Shen12_Cor_ad"]] = t2

t1 = read_tsv(file = paste0(in_dir,"/pombo_NPC_NcoI.enh_1e-2_mm10.bed"),col_names = FALSE, col_types = "cnncc") 
t2 = t1[,c(1:4)]
t2$X1 = paste0("chr",t2$X1)
t3 = as_tibble(str_split(t2$X4, ":", simplify = TRUE))
t2$X4 = t3$V1
feat_list[["Enh_Ron17_NPC"]] = t2

t1 = read_tsv(file = paste0(in_dir,"/mCO.enh.1e-2_mm10.bed"),col_names = FALSE, col_types = "cnncc") 
t2 = t1[,c(1:4)]
t2$X1 = paste0("chr",t2$X1)
t3 = as_tibble(str_split(t2$X4, ":", simplify = TRUE))
t2$X4 = t3$V1
feat_list[["Enh_Ron17_Cor_ad"]] = t2
##############


##############################################################################
# add annotations from reference datasets 
##############################################################################

### function: annotate_set

# identify overlap of elements from different peak and annotation files 

# input: peak_set (peaks with coordinates): columns: peak_id, Chr, Start, End, [...]
#        annot_set (coordinates of annotated elements): columns: Chr, Start, End, Feature (don't need to be labelled)

# output: peak_set + Feature from annot_set (if multiple annotations separated by "/")

annotate_set = function(peak_set, annot_set, chrom_set, annot_name){
  
  # annot_feature: returns features overlapping with individual region (define by vector with start and end coordinates)
  # features need to be filtered for relevant chromosome
  annot_feature = function(range_vector, feat_chrom){
    feat_elem = feat_chrom$Feature[(feat_chrom$Start >= range_vector[1] & feat_chrom$Start <= range_vector[2]) | 
                                     (feat_chrom$End >= range_vector[1] & feat_chrom$End <= range_vector[2])|
                                     (feat_chrom$Start <= range_vector[1] & feat_chrom$End >= range_vector[2])]
    new_Feature = paste0(feat_elem, collapse = "/")
    return(new_Feature)
  }
  
  #annotate_chrom: returns all annotations for specified chromosome; 
  #  filters peak_set and annot_set accordingly (annot_name defines column name for new annotation)
  annotate_chrom = function(Chr, peak_set, annot_set, annot_name){
    
    peaks_chrom = peak_set[peak_set$Chr == Chr,]
    feat_chrom = annot_set[annot_set$Chr == Chr,]
    
    new_Feature = apply(cbind(peaks_chrom$Start, peaks_chrom$End), 1, annot_feature, feat_chrom)
    peaks_chrom$new_Feature = new_Feature
    names(peaks_chrom)[names(peaks_chrom) == "new_Feature"] = annot_name
    
    return(peaks_chrom)
  }
  
  names(annot_set) = c("Chr", "Start", "End", "Feature")
  
  peak_set_list = lapply(chrom_set, annotate_chrom, 
                         peak_set = peak_set, annot_set = annot_set, annot_name = annot_name)
  
  peak_set_with_annot = bind_rows(peak_set_list)
  
  return(peak_set_with_annot)
}


### annotate features from all feature sets (in feat_list)
t1 = peak_annot1

for (i in 1:length(feat_list)){

  t1 = annotate_set(peak_set = t1, annot_set = feat_list[[i]], 
                    chrom_set = Chr, annot_name = names(feat_list)[i])
  message(paste0("Set ", i, " of ", length(feat_list), " annotated"))
}

peaks_full_annot = t1

write_csv(peaks_full_annot, path= paste0(out_dir,"/01_ATAC_peaks_ext_annot.csv"))

#########


############################################################################
# final peak classification
############################################################################
# default: "other"
# promoter assignment: all peaks assigned as promoter in original annotation or in Ensembl annotation 
# enhancers: Ensembl "enhancer" annotation or overlapping with H3K4m1 from reference datasets (and not annotated as promoter)

peaks_full_annot = read_csv(paste0(out_dir,"/01_ATAC_peaks_ext_annot.csv"))

t1 = as_tibble(peaks_full_annot)
t1$peak_type = "Other"
t1$linked_gene = t1$`Nearest PromoterID`

t1$peak_type[grepl("Enhancer", t1$Ensemble_feature_type)| 
               !is.na(t1$H3K4m1_Fo_P0)|!is.na(t1$H3K4m1_Cor_8w)] = "Enhancer"

t1$peak_type[grepl("promoter-TSS", t1$Annotation)] = "Promoter"

t3 = str_split(t1$Ensemble_feature_type,"/")
v1 = lapply(t3, `==`, "Promoter")
v2 = as_vector(lapply(v1, any))
t1$peak_type[v2] = "Promoter"

peak_types_annot = t1
#####

############################################################################
# distal element-promoter assignment 
#
#linked_gene: peaks assigned to closest TSS + regulatory regions from reference sets included for non-promoter elements
#linked_gene_100kb: peaks only assigned to closest TSS (max 100 kb cutoff)
############################################################################

gene_regulatory_region_columns = c(
  "Enh_Shen12_brain_E14","Enh_Shen12_Cor_ad",   
  "Enh_Ron17_NPC"  ,       "Enh_Ron17_Cor_ad" 
)

t2 = peak_types_annot 

for (i in gene_regulatory_region_columns){
  
  #all enhancers and "Other" elements with entry in regulatory column
  v1 = (t2$peak_type %in% c("Enhancer", "Other")) & !is.na(as_vector(t2[,i]))
  t3 = t2[v1,]
  t3$linked_gene = paste(t3$linked_gene,as_vector(t3[,i]), sep="/")
  t2[v1,] = t3
  
}

#remove duplicated gene names from final assignment
t3 = t2$linked_gene
t4 = str_split(t3,"/")
t5 = lapply(t4, unique)
t6 = lapply(t5, paste, collapse = "/")
t2$linked_gene = as_vector(t6)

#add bed peak_label: peak_id + peak type + linked genes
t2$peak_label = paste(t2$Peak_id, t2$peak_type, t2$linked_gene, sep = "_")

#linked_gene_100kb: peaks only assigned to closest TSS (max 100 kb cutoff)
t2$linked_gene_100kb[t2$`Distance to TSS`<100000] = t2$`Nearest PromoterID`[t2$`Distance to TSS`<100000]

# save annotation file and bed file
peaks_full_annot = t2
write_csv(peaks_full_annot, path= paste0(out_dir,"/01_ATAC_peaks_full_annot.csv"))

annot_bed = peaks_full_annot[,c("Chr", "Start", "End", "peak_label")]
write_tsv(annot_bed, path= paste0(out_dir,"/01_Peaks_impr_annot.bed"), col_names = FALSE)

########################


##############################################################################
##############################################################################
# [2] DESeq2 differential accessibility analysis and basic peak characterisation
##############################################################################
##############################################################################

##############################################################################
# load input (DESeq2 results from pipeline, group specifications and table and bed file with extended annotations)
##############################################################################

deseq_results0 = read_tsv(paste0(in_dir,"/merged_peaks.results.txt"))

gr_tab = read_csv(paste0(in_dir,"/group_table.csv"))
comp_tab = read_csv(paste0(in_dir,"/comparisons.csv"))

atac_bed0 = read_tsv(paste0(out_dir,"/01_Peaks_impr_annot.bed"), col_names = FALSE)
atac_annot = read_csv(paste0(out_dir,"/01_ATAC_peaks_full_annot.csv"))

#define chromosome set (to remove peaks in uncompletely annotated regions)
Chr = paste0("chr", c(1:19, "X", "Y"))

#define cutoff normalised count (.pseudo), log2FC, p

min_pseudo = 25
min_log2FC = 1
max_p = 0.05
##########

##############################################################################
# Run DESeq2
##############################################################################


coldata1 = as.data.frame(cbind(group = gr_tab$group), row.names = gr_tab$sample)

#get count matrix from original DESeq dataset, re-order according to sample order in group_tab0
t1 = deseq_results0
cts1 = t1[,match(paste0(gr_tab$sample, ".raw"), colnames(t1))]
names(cts1) = gr_tab$sample

#run DESeq
dds = DESeqDataSetFromMatrix(countData = cts1,
                             colData = coldata1,
                             design = ~ group)
dds = DESeq(dds)

save(dds, file = paste0(out_dir,"/02_Deseq_dataset.rda"))
############

##############################################################################
# Extract DESeq2 data and clean peak set, identify differentially accessible peaks
##############################################################################

#extract raw counts and normalized counts (.pseudo)

raw.counts = counts(dds,normalized=FALSE)
colnames(raw.counts) = paste(colnames(raw.counts),'raw',sep='.')

pseudo.counts = counts(dds,normalized=TRUE)
colnames(pseudo.counts) = paste(colnames(pseudo.counts),'pseudo',sep='.')


# get log2-transformed expression data 

ntd = normTransform(dds)
log_matrix =  as_tibble(assay(ntd))
names(log_matrix) = paste(colnames(log_matrix),'log',sep='.')

expr_tab = as_tibble(cbind(deseq_results0$Geneid, raw.counts,pseudo.counts, log_matrix))
names(expr_tab)[1] = "Peak_id"


#extract statistics for specified comparisons (in comparisons.csv), add to expr_tab

deseq_results_tab = expr_tab
i=1
for (i in 1:nrow(comp_tab)){
  t2 = results(dds, contrast = c("group", comp_tab$group1[i], 
                                 comp_tab$group2[i]))
  t3 = as_tibble(t2@listData)
  names(t3) = paste(comp_tab$group1[i],"vs",comp_tab$group2[i],names(t3),sep=".")
  deseq_results_tab = cbind(deseq_results_tab, t3)
}

deseq_results_tab = deseq_results_tab[order(deseq_results_tab$Peak_id),]


#add peak annotations, save
t1 = atac_annot[order(atac_annot$Peak_id),]

deseq_annot = as_tibble(cbind(t1[,c("Peak_id", "Chr", "Start", "End", "Strand", 
                                    "peak_type", "linked_gene", "peak_label", "linked_gene_100kb")],
                              deseq_results_tab[match(t1$Peak_id, deseq_results_tab$Peak_id),-1]))

write_csv(deseq_annot, file = paste0(out_dir,"/02_Deseq_results.csv"))


### remove low count peaks (normalized count (.pseudo) below cutoff (25) in all samples), and peaks not annotated on canonical chromosomes
t1 = deseq_annot              
gr_cols = paste0(gr_tab$sample,".pseudo")
pseudo_max = apply(t1[,paste0(gr_tab$sample,".pseudo")],1,max)

t2 = t1[ (t1$Chr %in% Chr) & (pseudo_max >= min_pseudo),]

deseq_cleaned = t2


#filter peaks with padj <0.05, abs(log2FC) > 1 in any pairwise group comparison 

t1 = deseq_cleaned
comparisons_padj_names = paste(comp_tab$group1,"vs",comp_tab$group2,"padj",sep=".")
comparisons_log2FC_names = paste(comp_tab$group1,"vs",comp_tab$group2,"log2FoldChange",sep=".")
comp_sign_diff = as_tibble((t1[,comparisons_padj_names] < max_p)&(abs(t1[,comparisons_log2FC_names]) > min_log2FC))
t2 = t1[apply(comp_sign_diff,1,any)&!apply(comp_sign_diff,1,anyNA),]
t3 = t1[apply(comp_sign_diff,1,anyNA),]

deseq_diff = t2

#generate corresponding bed files
t1 = atac_bed0

bed_cleaned = t1[t1$X4 %in% deseq_cleaned$peak_label,]
bed_diff = t1[t1$X4 %in% deseq_diff$peak_label,]

#save deseq and bed files for annotated/differential/gained/lost peaks

write_csv(deseq_annot, file = paste0(out_dir,"/02_deseq_annot.csv"))
write_csv(deseq_cleaned, file = paste0(out_dir,"/02_deseq_cleaned.csv"))
write_csv(deseq_diff, file = paste0(out_dir,"/02_deseq_diff.csv"))

write_delim(bed_cleaned, file = paste0(out_dir,"/02_peaks_cleaned.bed"), delim = "\t", col_names = FALSE)
write_delim(bed_diff, file = paste0(out_dir,"/02_peaks_diff.bed"), delim = "\t", col_names = FALSE)

#generate peak_type stat
t1 = deseq_annot %>% group_by(peak_type) %>% summarize(N_all_peaks = n())
t2 = deseq_cleaned %>% group_by(peak_type) %>% summarize(N_cleaned_peaks = n())
t3 = deseq_diff %>% group_by(peak_type) %>% summarize(N_diff_peaks = n())

peak_stat = cbind(t1, N_cleaned_peaks = t2$N_cleaned_peaks, N_diff_peaks = t3$N_diff_peaks)
write_csv(peak_stat, file = paste0(out_dir,"/02_peak_type_stats.csv"))

#############

##############################################################################
# additional comparisons: regions opening in mature astrocytes (Astr_P4 => Astr_2m), 
# lacking accessibility in vitro (Astr_BMP vs Astr_2m)
##############################################################################
deseq_diff = read_csv(file = paste0(out_dir,"/02_deseq_diff.csv"))
bed_diff = read_delim(file = paste0(out_dir,"/02_peaks_diff.bed"), delim = "\t", col_names = FALSE)

t1 = deseq_diff

t2 = t1[t1$Astr_P4.vs.Astr_2m.log2FoldChange < -1 & t1$Astr_P4.vs.Astr_2m.padj < 0.05 &
          t1$Astr_BMP.vs.Astr_2m.log2FoldChange < -1 & t1$Astr_BMP.vs.Astr_2m.padj < 0.05,]

t3 = bed_diff[bed_diff$X4 %in% t2$peak_label,]

write_csv(t2, file = paste0(out_dir,"/02_deseq_mat_opening_BMP_closed.csv"))
write_delim(t3, file = paste0(out_dir,"/02_peaks_mat_opening_BMP_closed.bed"), delim = "\t", col_names = FALSE)

#############




##############################################################################
##############################################################################
# [3] Integration with RNA-Seq data
##############################################################################
##############################################################################

##############################################################################
# load RNA and ATAC datasets, group tables and comparison tables
##############################################################################

#load raw input: bed with ATAC peaks, peak annotations, DESeq results, experimental design file
atac_bed0 = read_tsv(paste0(out_dir,"/02_peaks_cleaned.bed"), col_names = FALSE)

atac_results0 = read_csv(paste0(out_dir,"/02_deseq_cleaned.csv"))
rna_results0 = read_csv(paste0(in_dir,"/01_Deseq_results_expressed.csv"))

atac_gr_tab0 = read_csv(paste0(in_dir,"/group_table.csv"))
rna_gr_tab0 = read_csv(paste0(in_dir,"/group_tab_RNA.csv"))

atac_comp_tab0 = read_csv(paste0(in_dir,"/comparisons.csv"))
atac_comp0 = paste0(atac_comp_tab0$group1,".vs.", atac_comp_tab0$group2)

rna_comp_tab0 = read_csv(paste0(in_dir,"/comparisons_RNA.csv"))
rna_comp0 = paste0(rna_comp_tab0$group1,".vs.", rna_comp_tab0$group2)
#########

##############################################################################
# integrate RNA and ATAC datasets
##############################################################################

# extract from ATAC results: peak_id, linked gene, 
#         relevant stats (for each comparison ".padj", ".log2FoldChange", 
#         replace by ".atac_padj", ".atac_log2FC"),
#   pseudoexpression for all samples (replace ".pseudo" by ".atac")

t1 = atac_results0
t2.1 = t1[,c("Peak_id", "linked_gene", "linked_gene_100kb", "peak_type", "peak_label")]
t2.2 = t1[,c(paste0(atac_comp0,".padj"), paste0(atac_comp0,".log2FoldChange"))]
names(t2.2) = c(paste0(atac_comp0,".atac_padj"), paste0(atac_comp0,".atac_log2FC"))
t2.3 = t1[,paste0(atac_gr_tab0$sample,".pseudo")]
names(t2.3) = paste0(atac_gr_tab0$sample,".atac")

atac_results = cbind(t2.1, t2.2, t2.3)

# extract from RNA results: gene, 
#   relevant stats (for each comparison ".padj", ".log2FoldChange", 
#                   replace by ".rna_padj", ".rna_log2FC"),
#   pseudoexpression for all samples (replace ".pseudo" by ".rna")

t1 = rna_results0
t2.1 = t1[,"gene"]
t2.2 = t1[,c(paste0(rna_comp0,".padj"), paste0(rna_comp0,".log2FoldChange"))]
names(t2.2) = c(paste0(rna_comp0,".rna_padj"), paste0(rna_comp0,".rna_log2FC"))
t2.3 = t1[,paste0(rna_gr_tab0$sample,".pseudo")]
names(t2.3) = paste0(rna_gr_tab0$sample,".rna")

rna_results = cbind(t2.1, t2.2, t2.3)

### split atac dataset to single row for each linked gene

t1 = atac_results

linked_genes = str_split(t1$linked_gene,"/")
names(linked_genes) = t1$Peak_id

v3 = rep(names(linked_genes), lengths(linked_genes))
t3 = as_tibble(cbind(Peak_id = v3, linked_gene = unlist(linked_genes)))

t4 = t1[match(t3$Peak_id, t1$Peak_id),]
t4$linked_gene = t3$linked_gene

atac_results_by_gene = t4


### merge RNA/ATAC tables

r1 = rna_results
a1 = atac_results_by_gene

t1 = cbind(a1, r1[match(a1$linked_gene, r1$gene),])
r2 = r1[!(r1$gene %in% a1$linked_gene),]

t2 = as_tibble(matrix(ncol = ncol(t1), nrow = nrow(r2)) ) 
names(t2) = names(t1)

t2[,names(r2)] = r2

t3 = rbind(t1, t2)

atac_rna_results_combined = t3[order(t3$gene),]

write_csv(atac_rna_results_combined, path = paste0(out_dir,"/03_ATAC_RNA_results_merged.csv"))

### remove not expressed genes and minor peaks
#   threshold RNA >= 10 norm counts in any sample or ATAC >= 25

t1 = atac_rna_results_combined

rna_count_cols = paste0(rna_gr_tab0$sample,".rna")
atac_count_cols = paste0(atac_gr_tab0$sample,".atac")

max_rna = apply(t1[,rna_count_cols],1,max)
max_atac = apply(t1[,atac_count_cols],1,max)

t1[(!is.na(max_atac) & max_atac < 25), names(atac_results)] = NA

t2 = t1[ (max_rna >=10 & !is.na(max_rna)) | 
           (max_atac >=25 & !is.na(max_atac)) ,]

atac_rna_results_detected = t2

write_csv(atac_rna_results_detected, path = paste0(out_dir,"/03_ATAC_RNA_merged_detected.csv"))
##########

##############################################################################
##############################################################################
# [4] analyse chromatin accessibility changes at all peaks linked to immature/mature genes 
##############################################################################
##############################################################################

##############################################################################
# manually define selected comparisons
##############################################################################

#load combined ATAC/RNA dataset and gene sets to analyse (immature/mature genes)

atac_rna_results_detected = read_csv(file = paste0(out_dir,"/03_ATAC_RNA_merged_detected.csv"))
GOI_tab = read_csv(paste0(in_dir,"/07_Cor vs Str immat mat genes_genes.csv"))


##### manually define selected comparisons
#extract comparisons 
#transform changes of group comparisons to binary values (up/not changed/down = -1/0/1)
#exclude genes without linked peaks and peaks without linked genes

t1 = atac_rna_results_detected

v1 = names(t1)[grepl(".atac_padj", names(t1))]
v2 = str_remove(v1, ".atac_padj")
comp = v2

atac_comp = paste0("atac_reg_", comp)
rna_comp = paste0("rna_reg_", comp)

t2 = as_tibble(matrix(nrow = nrow(t1), ncol = 2*length(comp)))
names(t2) = c(atac_comp, rna_comp)
t2[,] = 0

for(i in 1:length(comp)){
  pa = t1[,paste0(comp[i],".atac_padj")]
  la = t1[,paste0(comp[i],".atac_log2FC")]
  t2[ (( pa <= 0.05 & la <= -1) & !(is.na(pa) | is.na(la)) ) , atac_comp[i]] = -1
  t2[ (( pa <= 0.05 & la >= 1) & !(is.na(pa) | is.na(la)) ) , atac_comp[i]] = 1
  pa = t1[,paste0(comp[i],".rna_padj")]
  la = t1[,paste0(comp[i],".rna_log2FC")]
  t2[ (( pa <= 0.05 & la <= -1) & !(is.na(pa) | is.na(la)) ) , rna_comp[i]] = -1
  t2[ (( pa <= 0.05 & la >= 1) & !(is.na(pa) | is.na(la)) ) , rna_comp[i]] = 1
}

t3 = as_tibble(cbind(t1,t2))

t4 = t3[!(is.na(t3$Peak_id) | is.na(t3$gene)),]

res_bin = t4

write_csv(res_bin, path = paste0(out_dir,"/03_ATAC_RNA_merged_detected_with_binary_comparisons.csv"))

# define groups

res_list = NULL

t1 = res_bin
t2 = t1[ t1$gene %in% na.omit(GOI_tab$mat_common),]
t3 = t2[t2$rna_reg_Astr_BMP.vs.Astr_2m <0 , ]
res_list[["mature_genes_low"]] = t3
res_list[["mature_genes_low_mature_peaks_closed"]] = t3[t3$atac_reg_Astr_BMP.vs.Astr_2m<0,]
t3 = t2[t2$rna_reg_Astr_BMP.vs.Astr_2m >=0 , ]
res_list[["mature_genes_high"]] = t3
res_list[["mature_genes_high_mature_peaks_closed"]] = t3[t3$atac_reg_Astr_BMP.vs.Astr_2m<0,]


#create statistics of number of genes per group

t1 = NULL

for (i in names(res_list)){
  t2 = tibble(group = i, N_genes = length(unique(na.omit(res_list[[i]]$gene))))
  t1 = rbind(t1, t2)
}

write_csv(t1, file = paste0(out_dir,"/04_stats_immat mat genes with vs without diff peaks.csv"))


#########


##############################################################################
# count total number of peaks lacking accessibility per gene
##############################################################################

t1 = res_bin

#determine number of total/opening/closing peaks by gene
t2 = t1 %>% group_by(gene) %>% summarise(N_peaks=n(), 
                                         N_mat_open_BMP_closed = 
                                           length(atac_reg_Astr_P4.vs.Astr_2m[
                                             atac_reg_Astr_BMP.vs.Astr_2m == -1]),
                                         N_mat_opening_BMP_closed = 
                                           length(atac_reg_Astr_P4.vs.Astr_2m[
                                             atac_reg_Astr_P4.vs.Astr_2m == -1 & 
                                               atac_reg_Astr_BMP.vs.Astr_2m == -1]))

peaks_reg_per_gene = t2

write_csv(peaks_reg_per_gene, file = paste0(out_dir,"/04_peaks_reg_per_gene_BMP_closed.csv"))

##############################################################################
# analyse and plot total number of peaks lacking accessibilty per gene in mature genes
##############################################################################

# frequency distribution plot table

t1 = peaks_reg_per_gene

t2 = t1[t1$gene %in% res_list$mature_genes_low$gene, ]
t2$group = "mature_genes_low"

t3 = t1[t1$gene %in% res_list$mature_genes_high$gene, ]
t3$group = "mature_genes_high"

pl_tab = rbind(t2,t3)

write_csv(pl_tab, file = paste0(out_dir,"/04_peaks_BMP_closed_mat_genes_BMP_low_vs_high_genes.csv"))


### plot peaks too low per gene

g1 = ggplot(data = pl_tab) + geom_boxplot(aes(x=group, y = N_mat_open_BMP_closed, fill = group),
                                      outlier.size = 0.5) +
  theme_light()+
  labs(x = "", y = "peaks_per_gene", fill = "")+
  theme(text =  element_text(size = 12, face = 2),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(file = paste0(out_dir,"/04_peaks_BMP_closed_mat_genes_BMP_low_vs_high_genes.pdf"), width = 3.2, height = 3)
{
  g1
}    
dev.off()

sink(paste0(out_dir,"/04_peaks_BMP_closed_mat_genes_BMP_low_vs_high_genes_wilcoxon_stat.txt"))
{
  pl_tab %>% group_by(group) %>% summarise(n())
  pairwise.wilcox.test(pl_tab$N_mat_open_BMP_closed, pl_tab$group,
                       p.adjust.method = "BH")
}
sink()

##############


