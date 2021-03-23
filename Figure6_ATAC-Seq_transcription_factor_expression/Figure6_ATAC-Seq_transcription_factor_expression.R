##############################################################################
# global settings
##############################################################################

in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(tidyverse)
library(DESeq2)
library(heatmaps)
library(colorRamps)
library(GenomicRanges)
library(rtracklayer)

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
##############################################################################
# [3] Downstream analysis
#
#
# define group comparisons
# extract coverage matrices from bigwig files and plot coverage heatmaps
##############################################################################
##############################################################################

deseq_cleaned = read_csv(file = paste0(out_dir,"/02_deseq_cleaned.csv"))
bed_cleaned = read_delim(file = paste0(out_dir,"/02_peaks_cleaned.bed"), delim = "\t", col_names = FALSE)

gr_tab = read_csv(paste0(in_dir,"/group_table.csv"))
bigwig_tab = read_csv(paste0(in_dir,"/group_table_bigwigs.csv"))

comp_tab = read_csv(paste0(in_dir,"/comparisons.csv"))


##############################################################################
# define group comparisons
##############################################################################

comp = paste0(comp_tab$group1,".vs.", comp_tab$group2)

t1 = deseq_cleaned

### create binarised accesibility table for each comparison -1/0/1 for increased/unchanged/decreased accessibility

m1 = matrix(data = 0, ncol = length(comp), nrow = nrow(t1), 
            dimnames = list(t1$Peak_id, comp))

for (i in comp){
  m1[(t1[paste0(i,".padj")] <= 0.05) & (t1[paste0(i,".log2FoldChange")] >= 1)  ,i] = 1
  m1[(t1[paste0(i,".padj")] <= 0.05) & (t1[paste0(i,".log2FoldChange")] <= -1)  ,i] = -1
}

cb = as_tibble(cbind(Peak_id = rownames(m1), m1))


### define peak modules, add Peak IDs to list

l1 = list()

#all peaks opening/closing upon TF OE
l1[["Rorb_opening"]] = cb$Peak_id[cb$EGFP.vs.Rorb == -1]
l1[["Rorb_closing"]] = cb$Peak_id[cb$EGFP.vs.Rorb == 1]
l1[["Dbx2_opening"]] = cb$Peak_id[cb$EGFP.vs.Dbx2 == -1]
l1[["Dbx2_closing"]] = cb$Peak_id[cb$EGFP.vs.Dbx2 == 1]
l1[["Lhx2_opening"]] = cb$Peak_id[cb$EGFP.vs.Lhx2 == -1]
l1[["Lhx2_closing"]] = cb$Peak_id[cb$EGFP.vs.Lhx2 == 1]

#only mature opening peaks
v1 = cb$Peak_id[cb$Astr_P4.vs.Astr_2m == -1]
l2 = lapply(l1, intersect, v1)
names(l2) = paste0("mature_opening_", names(l1))
l2[["mature_opening"]] = v1

#only mature closing peaks
v1 = cb$Peak_id[cb$Astr_P4.vs.Astr_2m == 1]
l3 = lapply(l1, intersect, v1)
names(l3) = paste0("mature_closing_", names(l1))
l3[["mature_closing"]] = v1

peak_list = c(l1,l2,l3)
##############

##############################################################################
# save module stats and module bed files
##############################################################################

### peak module stats

peak_module_stats = tibble(module = names(peak_list), N_peaks = lengths(peak_list))
write_csv(peak_module_stats, file = paste0(out_dir,"/03_peak_module_stats.csv"))

### save bed files for each peak module

module_bed_dir = paste0(out_dir,"/03_module_bed_files/")

if (!dir.exists(module_bed_dir)) {dir.create(module_bed_dir)}

module_bed_list = peak_list

t1 = bed_cleaned

for (i in names(peak_list)){
  
  peak_label =  deseq_cleaned$peak_label[deseq_cleaned$Peak_id %in% peak_list[[i]]]
  module_bed_list[[i]] = t1[t1$X4 %in% peak_label,]
  
  write_tsv(module_bed_list[[i]], file = paste0(module_bed_dir,"/03_Peaks_",i,".bed"), col_names = FALSE)
  
}


#############


##############################################################################
# create and plot coverage matrices for selected comparisons
##############################################################################

bigwig_group_tab = read_csv(paste0(in_dir,"/group_table_bigwigs.csv"))
peak_coord = import(paste0(out_dir,"/02_peaks_cleaned.bed"))

#convert peak IDs to peak labels of bed files
peak_label_list = peak_list

for (i in names(peak_label_list)){
  peak_label_list[[i]] = deseq_cleaned$peak_label[deseq_cleaned$Peak_id %in% peak_list[[i]]]
}

#define peak sets to plot
names(peak_label_list)
peak_label_list_plot = peak_label_list[c("Rorb_opening", "Lhx2_opening")]

#get coordinates for peak center +/- 1000 bp
peak_coord_ext = peak_coord
peak_center = (start(peak_coord)+end(peak_coord))%/%2
start(peak_coord_ext) = peak_center - 1000
end(peak_coord_ext) = peak_center + 1000


### for each bigwig file in bigwig_group_tab, extract coverage matrix for each specified peak set 
#   (each 100 bins, +/- 1000 bp from peak center)
#   takes some time! (has to load and go through all bigwig files)

#function: create frip-normalised covarage heatmap for selected peaks fo peak_coord set from bigwig file

get_norm_cov_heatmap = function(peak_labels, peak_coord, bigwig, frip = 1){
  
  #extract peak coordinates for specific set
  t2 = peak_coord[peak_coord$name %in% peak_labels,]
  
  #create coverage heatmap object and normalise by fraction in peaks (frip) score for specific sample
  t3 = CoverageHeatmap(windows = t2, track = bigwig1,
                       coords = c(-1000, 1000), label = "", nbin = 100)
  t3@image/frip
  t3@scale/frip
  return(t3)
}


# get heatmap matrices for each peakset in each bigwig file 
#     (pl_matrix_list = list[bigwig_group_tab$group] >> list[peak_label_list_plot] >> Heatmap object for peaks in peak_label_list_plot)

cov_matrix_list = list()

for (i in 1:nrow(bigwig_group_tab)){
  
  message(paste0("calculating matrices for sample ", i, " of ",nrow(bigwig_group_tab)))
  
  bigwig1 = import(paste0(bigwig_group_tab$path_bigwig[i],"/", bigwig_group_tab$file_bigwig[i]),
                   format="BigWig", as = "RleList")
  
  cov_matrix_list[[bigwig_group_tab$group[i]]] = lapply(peak_label_list_plot, get_norm_cov_heatmap, 
                                                       peak_coord = peak_coord_ext, 
                                                       bigwig = bigwig1, 
                                                       frip = bigwig_group_tab$frip[i])
}

save(cov_matrix_list, peak_label_list_plot, bigwig_group_tab, 
     file = paste0(out_dir,"/03_coverage_matrices.RData"))
#############


####################################################################################
# plot coverage matrices
####################################################################################

load(file = paste0(out_dir,"/03_coverage_matrices.RData"))

### for plotting: 
# sort all plot matrices into list of list of matrices to be plotted side by side
# (log2+1) transform coverage values and adjust scale and order to match reference matrix (Astr_2m)
# plot matrices from one peak module for all groups combined

pl_matrix_list = list()
ref_matrix_name = "Astr_2m"

for (i in names(peak_label_list_plot)){
  
  l1 = list()
  
  ref_mat = cov_matrix_list[[ref_matrix_name]][[i]]
  t1 = log2(ref_mat@image + 1)
  ref_sum_int = apply(t1, 1, sum)
  ref_scale = c(0, 0.8*max(t1))
  
  for (j in bigwig_group_tab$group){
    m1 = cov_matrix_list[[j]][[i]]
    m1@image = log2(m1@image + 1)
    m1@image = m1@image[order(-ref_sum_int),]
    m1@scale = ref_scale
    pl_matrix_list[[i]][[j]] = m1 
  }
  
  pdf(file = paste0(out_dir,"/03_coverage_heatmaps_",i,".pdf"), 
      width = length(pl_matrix_list[[i]])*5, height = length(ref_sum_int)/40 +0.5)
  {
    plotHeatmapList(pl_matrix_list[[i]], color = matlab.like(50), legend = FALSE)
    plotHeatmapList(pl_matrix_list[[i]], color = matlab.like(50), legend = TRUE)
    
  }
  dev.off()
  
}




