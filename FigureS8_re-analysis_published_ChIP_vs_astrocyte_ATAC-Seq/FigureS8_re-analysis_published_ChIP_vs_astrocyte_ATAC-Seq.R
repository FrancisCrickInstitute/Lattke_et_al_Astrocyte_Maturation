##############################################################################
# global settings
##############################################################################

in_dir = "./00_input"
out_dir = "./01_output"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(tidyverse)

###############


##############################################################################
##############################################################################
# [1] peak annotations from ATAC dataset with reference datasets 
##############################################################################
##############################################################################


##############################################################################
# load ATAC peak annotations from pipeline and reference datasets for annotations 
##############################################################################

peak_annot1 = read_csv(paste0(in_dir,"/01_ATAC_peaks_full_annot.csv"))

Chr = paste0("chr",c(1:19, "X", "Y"))

### load files with features for additional annotations

#extract location and peak name from Fezf2 ChIP pipeline "[replicate]_peaks.annotatePeaks.txt"

extract_chip_peaks = function(peak_file){
  t1 = read_tsv(file = peak_file)
  t2 = t1[,c(2,3,4,1)]
  names(t2) = c("Chr", "Start", "End", "Feature") 
  t2$Chr = paste0("chr", t2$Chr)
  return(t2)
}

feat_list = list()

feat_list[["ATAC_mature_opening"]] = read_tsv(file = paste0(in_dir,"/03_Peaks_mature_opening.bed"), col_names = FALSE)
feat_list[["ATAC_Rorb_opening"]] = read_tsv(file = paste0(in_dir,"/03_Peaks_Rorb_opening.bed"), col_names = FALSE)
feat_list[["ATAC_Lhx2_opening"]] = read_tsv(file = paste0(in_dir,"/03_Peaks_Lhx2_opening.bed"), col_names = FALSE)

feat_list[["Fezf2_R1_Lodato14"]] = extract_chip_peaks(paste0(in_dir,"/Fezf2_R1_peaks.annotatePeaks.txt"))
feat_list[["Fezf2_R2_Lodato14"]] = extract_chip_peaks(paste0(in_dir,"/Fezf2_R2_peaks.annotatePeaks.txt"))
feat_list[["RORg_R1_He17"]] = extract_chip_peaks(paste0(in_dir,"/RORgt_R1_peaks.annotatePeaks.txt"))
feat_list[["RORg_R2_He17"]] = extract_chip_peaks(paste0(in_dir,"/RORgt_R2_peaks.annotatePeaks.txt"))
feat_list[["Lhx2_R1_Monahan17"]] = extract_chip_peaks(paste0(in_dir,"/Lhx2_chip_R1_peaks.annotatePeaks.txt"))
feat_list[["Lhx2_R2_Monahan17"]] = extract_chip_peaks(paste0(in_dir,"/Lhx2_chip_R2_peaks.annotatePeaks.txt"))



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

t1[t1 == ""] = NA

peaks_full_annot = t1

write_csv(peaks_full_annot, file = paste0(out_dir,"/01_ATAC_peaks_ext_annot.csv"))

save(peaks_full_annot, feat_list, file = paste0(out_dir,"/01_ATAC_peaks_ext_annot_and_annot_sets.rda"))

#########



##############################################################################
##############################################################################
# [2] identify and quantify overlaps between ATAC dataset and reference peak sets
##############################################################################
##############################################################################

##############################################################################
# filter peaks overlapping with specific set into list and create overlap stats 
##############################################################################

load(file = paste0(out_dir,"/01_ATAC_peaks_ext_annot_and_annot_sets.rda"))

overlap_stat = tibble(comparison = names(feat_list), 
                      N_peaks_total = as_vector(lapply(feat_list, nrow)), 
                      N_overlap = 0)

overlap_peak_list = list()

for (i in names(feat_list)){
  overlap_peak_list[[i]] = peaks_full_annot$Peak_id[!is.na(as_vector(peaks_full_annot[,i]))]
  overlap_stat$N_overlap[overlap_stat$comparison == i] = length(overlap_peak_list[[i]])
}

write_csv(overlap_stat, file = paste0(out_dir,"/02_Stat_ATAC_peakset_overlap_with_reference_sets.csv"))
###########

##############################################################################
# manually add additional overlap comparisons
##############################################################################

l1 = overlap_peak_list

l1[["Fezf2_Lodato14_R1_R2_overlap"]] = intersect(l1$Fezf2_R1_Lodato14, l1$Fezf2_R2_Lodato14)
l1[["RORg_He17_R1_R2_overlap"]] = intersect(l1$RORg_R1_He17, l1$RORg_R2_He17)
l1[["Lhx2_Monahan17_R1_R2_overlap"]] = intersect(l1$Lhx2_R1_Monahan17, l1$Lhx2_R2_Monahan17)
l1[["RORg_He17_Fezf2_Lodato14_overlap"]] = intersect(l1$Fezf2_Lodato14_R1_R2_overlap, l1$RORg_He17_R1_R2_overlap)

l1[["ATAC_RORb_opening_vs_RORg_He17_R1_R2_overlap"]] = intersect(l1$ATAC_Rorb_opening, l1$RORg_He17_R1_R2_overlap )
l1[["ATAC_RORb_opening_vs_Fezf2_Lodato14_R1_R2_overlap"]] = intersect(l1$ATAC_Rorb_opening, l1$Fezf2_Lodato14_R1_R2_overlap)
l1[["ATAC_RORb_opening_vs_Lhx2_Monahan17_R1_R2_overlap"]] = intersect(l1$ATAC_Rorb_opening, l1$Lhx2_Monahan17_R1_R2_overlap)
l1[["ATAC_Lhx2_opening_vs_RORg_He17_R1_R2_overlap"]] = intersect(l1$ATAC_Lhx2_opening, l1$RORg_He17_R1_R2_overlap )
l1[["ATAC_Lhx2_opening_vs_Fezf2_Lodato14_R1_R2_overlap"]] = intersect(l1$ATAC_Lhx2_opening, l1$Fezf2_Lodato14_R1_R2_overlap)
l1[["ATAC_Lhx2_opening_vs_Lhx2_Monahan17_R1_R2_overlap"]] = intersect(l1$ATAC_Lhx2_opening, l1$Lhx2_Monahan17_R1_R2_overlap)
l1[["ATAC_mature_opening_vs_RORg_He17_R1_R2_overlap"]] = intersect(l1$ATAC_mature_opening, l1$RORg_He17_R1_R2_overlap)
l1[["ATAC_mature_opening_vs_Fezf2_Lodato14_R1_R2_overlap"]] = intersect(l1$ATAC_mature_opening, l1$Fezf2_Lodato14_R1_R2_overlap)
l1[["ATAC_mature_opening_vs_Lhx2_Monahan17_R1_R2_overlap"]] = intersect(l1$ATAC_mature_opening, l1$Lhx2_Monahan17_R1_R2_overlap)

overlap_stat_ext = tibble(comparison = names(l1), N_peaks = lengths(l1))
write_csv(overlap_stat_ext, file = paste0(out_dir,"/02_Stat_ATAC_peakset_overlap_with_reference_sets_extended.csv"))

overlap_peak_list_ext = l1

save(peaks_full_annot, overlap_peak_list, overlap_peak_list_ext,
     file = paste0(out_dir,"/02_peak_overlap_analysis.rda"))


#save bed files with Rorb-opening RORg binding and Lhx2-opening Lhx2 binding peaks

t1 = feat_list$ATAC_Rorb_opening
t2 = peaks_full_annot[peaks_full_annot$Peak_id %in% overlap_peak_list_ext$ATAC_RORb_opening_vs_RORg_He17_R1_R2_overlap,]
t3 = t1[t1$X4 %in% t2$ATAC_Rorb_opening,]

write_tsv(t3,file = paste0(out_dir,"/02_Rorb_opening_RORg_binding_peaks.bed"), col_names = FALSE)

t1 = feat_list$ATAC_Lhx2_opening
t2 = peaks_full_annot[peaks_full_annot$Peak_id %in% overlap_peak_list_ext$ATAC_Lhx2_opening_vs_Lhx2_Monahan17_R1_R2_overlap,]
t3 = t1[t1$X4 %in% t2$ATAC_Lhx2_opening,]

write_tsv(t3,file = paste0(out_dir,"/02_Lhx2_opening_Lhx2_binding_peaks.bed"), col_names = FALSE)



##########



##############################################################################
##############################################################################
# [3] integrated analysis with RNA-Seq data
##############################################################################
##############################################################################

load( file = paste0(out_dir,"/02_peak_overlap_analysis.rda"))

### load RNA-Seq datasets (atac_rna_ref_set: merged set with annotated peaks) 
# gene of interest sets (from combined TF OE experiment)

atac_rna_ref_set = read_csv(file = paste0(in_dir,"/03_ATAC_RNA_merged_detected_with_binary_comparisons.csv") )

GOI_tab = read_csv(file = paste0(in_dir,"/02_genes_in_modules.csv") )

# extract specific gene sets

l1 = lapply(as.list(GOI_tab), na.omit)
names(l1)
l2 = l1[c(
  "Rorb_up", "Fezf2_up", "comb_RF_not_R_F_up", "mature_low_in_EGFP_Rorb_up" ,               
  "mature_low_in_EGFP_Fezf2_up", "mature_low_in_EGFP_comb_RF_not_R_F_up"        
)]

l2[["mature_low_in_EGFP_not_induced"]] = l1$mature_low_in_EGFP[!(l1$mature_low_in_EGFP %in% unlist(l2) )]
l2[["mature_low_in_EGFP_all"]] = l1$mature_low_in_EGFP 

reg_genes_list = l2

# stats related gene sets

t1 = tibble(gene_set = names(reg_genes_list), N_genes = lengths(reg_genes_list), N_peaks_total =  0)
for (i in names(reg_genes_list)){
  t1$N_peaks_total[t1$gene_set == i] = 
    length(unique(atac_rna_ref_set$Peak_id[atac_rna_ref_set$gene %in% reg_genes_list[[i]] ]) )
}

write_csv(t1, file = paste0(out_dir,"/03_reg_gene_set_stats.csv"))


##############################################################################
# analyse specific gene - peak relationships
##############################################################################

### identify potential Fezf2 target genes
t1 = atac_rna_ref_set[,c("Peak_id","linked_gene", "peak_type", "peak_label")]
t2 = t1[t1$Peak_id %in% overlap_peak_list_ext$Fezf2_Lodato14_R1_R2_overlap,]

l1 = reg_genes_list

for (i in names(l1)){
  l1[[i]] = t2[t2$linked_gene %in% reg_genes_list[[i]],]
}

names(l1) = paste0("Fezf2_targets_", names(l1))

atac_peak_gene_comp_list = l1


### identify potential Rorb target genes

t1 = atac_rna_ref_set[,c("Peak_id","linked_gene", "peak_type", "peak_label")]
t2 = t1[t1$Peak_id %in% overlap_peak_list_ext$RORg_He17_R1_R2_overlap,]

l1 = reg_genes_list

for (i in names(l1)){
  l1[[i]] = t2[t2$linked_gene %in% reg_genes_list[[i]],]
}

names(l1) = paste0("RORg_targets_", names(l1))

atac_peak_gene_comp_list = c(atac_peak_gene_comp_list, l1)

### identify potential Rorb Fezf2 combined target genes

t1 = atac_rna_ref_set[,c("Peak_id","linked_gene", "peak_type", "peak_label")]
t2 = t1[t1$Peak_id %in% overlap_peak_list_ext$RORg_He17_Fezf2_Lodato14_overlap,]

l1 = reg_genes_list

for (i in names(l1)){
  l1[[i]] = t2[t2$linked_gene %in% reg_genes_list[[i]],]
}

names(l1) = paste0("RORg_Fezf2_cobind_targets_", names(l1))

atac_peak_gene_comp_list = c(atac_peak_gene_comp_list, l1)


###create gene number stats and peak-gene tables for comparisons

l1 = atac_peak_gene_comp_list

if (!dir.exists("./01_output/03_peak_gene_links_sel_comp/")) {dir.create("./01_output/03_peak_gene_links_sel_comp/")}

peak_gene_link_stats = NULL

for (i in names(l1)){
  t1 = l1[[i]]
  write_csv(t1, file = paste0(out_dir,"/03_peak_gene_links_sel_comp/",i,"_peak_gene_tab.csv"))
  t2 = tibble(comparison = i, N_peak_gene_links = nrow(t1), 
              N_peaks = length(unique(t1$Peak_id)), N_genes = length(unique(t1$linked_gene)))
  peak_gene_link_stats = rbind(peak_gene_link_stats, t2)
}

write_csv(peak_gene_link_stats, file = paste0(out_dir,"/03_peak_gene_links_stats.csv"))

##########





