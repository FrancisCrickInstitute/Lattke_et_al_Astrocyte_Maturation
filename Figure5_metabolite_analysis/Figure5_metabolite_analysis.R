library(tidyverse)
library(pheatmap)
library(colorRamps)

if (!dir.exists("./out/")) dir.create("./out/")

metab_tab = read_csv("metabolites_cells_and_media.csv")

metab_class_tab = read_csv("metabolite_classes.csv")

metab = metab_class_tab$Metabolite
exper = unique(metab_tab$Exp)
cond = c("EGFP", "Rorb", "Fezf2")


###############################################################
# calculate normalised data
###############################################################

t1 = metab_tab

#calculate concentration and labelling rel to mean of experiment (incl log2+1 transformed values)

t3 = NULL

for (i in metab){
  
  for (j in exper){
    t2 = t1[t1$Metabolite == i & t1$Exp == j,]
    t2$conc_norm_cells = t2$conc_cells/mean(t2$conc_cells)
    t2$conc_log_norm_cells = log2(t2$conc_cells+0.001) - mean(log2(t2$conc_cells+0.001))
    t2$conc_media_net = t2$conc_media - t2$conc_media_baseline
    t2$conc_log_norm_media = log2(t2$conc_media+0.001) - mean(log2(t2$conc_media+0.001))
    t3 = rbind(t3, t2)
  }
  
}

# create merged sample identifier
t3$sample = paste0(t3$Condition,"_", t3$Exp)

# add metabolite class
t3$Pathway_Category = metab_class_tab$Pathway_Category[match(t3$Metabolite, metab_class_tab$Metabolite)]

# reorder by condition, then metabolites 
t3 = t3[order(match(t3$Condition, cond)),]
t3 = t3[order(match(t3$Metabolite, metab)),]

metab_tab_ext = t3



#############################################
#define variables for plotting/analysis
#############################################

plot_var = c("conc_cells", "conc_media", "conc_media_net", "conc_log_norm_cells", 
             "conc_log_norm_media", "label_percent_cells", "label_percent_media")


#############################################
# statistics (series of t-tests for each metabolite vs EGFP control, BH correction for p values)
#############################################

t_test_vs_control = function(metab_tab_ext, var){
  
  t1 = metab_tab_ext[!is.na(as_vector(metab_tab_ext[, var])),]
  
  #pairwise t-tests for all metabolites EGFP vs Rorb/Fezf2 with global Benjamini-Hochberg correction
  t3 = NULL
  
  for (i in unique(t1$Metabolite)){
    t2 = t1[t1$Metabolite == i,]
    s2 = pairwise.t.test(as_vector(t2[, var]), t2$Condition, 
                         p.adjust.method = "none")
    s3 = tibble(Metabolite = i, 
                comp_EGFP_vs = rownames(s2$p.value), 
                p = as_vector(s2$p.value[,"EGFP"]))
    t3 = rbind(t3, s3)
  }
  
  t3$padj = p.adjust(t3$p, method = "BH")
  t3$sign = ""
  t3$sign[t3$padj <= 0.05] = "*"
  t3$sign[t3$padj <= 0.01] = "**"
  t3$sign[t3$padj <= 0.001] = "***"
  
  stat_tab = t3
  
  write_csv(stat_tab, file = paste0("./out/stat_", var,".csv"))
  
  return(stat_tab)
  
}

stat_list = lapply(plot_var, t_test_vs_control, metab_tab_ext = metab_tab_ext)
names(stat_list) =plot_var


######################################################
# individual dotplots for all metabolites (raw and normalised concentration and C-13 labelling fraction)
######################################################

overview_dotplot = function(metab_tab_ext, y, scales = "free"){
  
  t1 = metab_tab_ext
  names(t1)[names(t1) == y] = "y"
  
  p1 = ggplot(data = t1)+ 
    facet_wrap(~Metabolite, scales = scales)+
    geom_point(aes(x = Condition, y = y, color = Exp, shape = Exp)) + 
    scale_x_discrete(limits = cond) + labs(y = y) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  pdf(file = paste0("./out/dotplot_",y,".pdf"), width = 15, height = 20)
  {plot(p1)}
  dev.off()
  
}


plot_var1 = c("conc_cells", "conc_media", "conc_media_net")
lapply(plot_var1, overview_dotplot, metab_tab_ext = metab_tab_ext, scales = "free")

plot_var2 = c("conc_log_norm_cells", "conc_log_norm_media", 
              "label_percent_cells", "label_percent_media")
lapply(plot_var2, overview_dotplot, metab_tab_ext = metab_tab_ext, scales = "fixed")




####################################################################################
# plot as heatmaps, save plotted values as table sample vs metabolite
####################################################################################

### function to create heatmap plot matrix and save mtrix with stats

create_plot_matr_save_with_stats = function(plot_var, metab_tab_ext, stat_list = NULL){
  
  t1 = metab_tab_ext
  
  # rename plot variable to plot_var
  names(t1)[names(t1) == plot_var] = "plot_var"
  
  t1 = t1[!is.na(t1$plot_var),]
  
  #create plot matrix 
  
  m1 = matrix(nrow = length(unique(t1$Metabolite)), ncol = length(unique(t1$sample)),
              dimnames = list(unique(t1$Metabolite), unique(t1$sample)) )
  
  for (i in unique(t1$Metabolite)){
    m1[i,] = t1$plot_var[t1$Metabolite == i]
  }
  
  plotmatr = m1
  
  #save table with plot_var data by sample and metabolite (with stats if available)
  
  t2 = as_tibble(cbind(Metabolite = rownames(m1), m1))
  
  if (!is.null(stat_list)){
    if (plot_var %in% names(stat_list)){
      stat_tab = stat_list[[plot_var]]
      t2$padj_EGFP_vs_Rorb = stat_tab$padj[stat_tab$comp_EGFP_vs == "Rorb"]
      t2$padj_EGFP_vs_Fezf2 = stat_tab$padj[stat_tab$comp_EGFP_vs == "Fezf2"]
    }
  }
  
  write_csv(t2, file = paste0("./out/table_",plot_var,"_values_and_stats.csv"))
  
  return (plotmatr)
}


### heatmap function

expr_heatmap = function(matr, pl_genes = NULL, scale_range = NULL,
                        scale_colors = colorRampPalette(c("blue", "white", "red"))(250),
                        cellwidth = 16, cellheight = 16, 
                        pl_file = "heatmap.pdf",
                        pl_title = pl_file){
  
  if (nrow(matr)>=2){
    
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
               color = scale_colors,
               breaks = seq(scale_range[1], scale_range[2], length.out = length(scale_colors)+1),
               border_color = NA, fontsize = cellheight, fontface = 4,
               cellwidth = cellwidth, cellheight = cellheight,
               main =  paste0(pl_title, " ( ", nrow(matr), " genes)")
      )
      
      pheatmap(matr, cluster_rows = FALSE, show_rownames=TRUE,
               cluster_cols=FALSE, 
               color = scale_colors,
               breaks = seq(scale_range[1], scale_range[2], length.out = length(scale_colors)+1),
               border_color = NA, fontsize = cellheight, fontface = 4,
               cellwidth = cellwidth, cellheight = cellheight,
               main =  paste0(pl_title, " (", nrow(matr), " genes)")
      )
      
    }
    dev.off()
  }
}

### generate plot matrices for all plot_var variables

matr_list = lapply(plot_var, create_plot_matr_save_with_stats, metab_tab_ext = metab_tab_ext, 
                  stat_list = stat_list)
names(matr_list) = plot_var

#extract medium baseline levels and add to 'table_conc_media_values_and_stats.csv' 

t1 = metab_tab_ext
t2 =t1[!duplicated(t1$Metabolite),]
t2 = t2[!is.na(t2$conc_media_baseline),]

t3 = read_csv("./out/table_conc_media_values_and_stats.csv")

#omit Tryptophan, no baseline level recorded
t3$conc_media_baseline = t2$conc_media_baseline[match(t3$Metabolite, t2$Metabolite)]
write_csv(t3, file = "./out/table_conc_media_values_and_stats.csv")


### plot heatmaps

# heatmaps for conc in cells (all/only signifcantly changed)

expr_heatmap(matr = matr_list$conc_log_norm_cells, pl_genes = NULL, scale_range = NULL,
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_conc_log_norm_cells.pdf")

v1 = unique(stat_list$conc_log_norm_cells$Metabolite[
  stat_list$conc_log_norm_cells$padj <= 0.05])

expr_heatmap(matr_list$conc_log_norm_cells, pl_genes = v1, scale_range = NULL,
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_conc_log_norm_cells_only_sign.pdf")


# heatmaps for conc in media (all/only signifcantly changed)

expr_heatmap(matr_list$conc_log_norm_media, pl_genes = NULL, scale_range = NULL,
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_conc_log_norm_media.pdf")

v1 = unique(stat_list$conc_log_norm_media$Metabolite[
  stat_list$conc_log_norm_media$padj <= 0.05])

expr_heatmap(matr_list$conc_log_norm_media, pl_genes = v1, scale_range = NULL,
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_conc_log_norm_media_only_sign.pdf")


# heatmaps for percent labelled in cells/media

expr_heatmap(matr_list$label_percent_cells, pl_genes = NULL, scale_range = c(0,100),
             scale_colors = matlab.like(250),
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_label_percent_cells.pdf")

expr_heatmap(matr_list$label_percent_media, pl_genes = NULL, scale_range = c(0,100),
             scale_colors = matlab.like(250),
             cellwidth = 16, cellheight = 16, 
             pl_file = "./out/heatmap_label_percent_media.pdf")





