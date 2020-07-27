rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas'
setwd(path)
# libraries
library(pacman)
p_load(pheatmap,grid,gridExtra,stringr,ggplot2, dplyr, tibble)

# function
plot_2_heatmap <- function(data, clean=FALSE, col_list, max_pro){
  # Packages required: pheatmap, grid, gridExtra, stringr
  if(missing(col_list)){
    col_count = apply(data,2,sum)
    col_count_sorted = sort(col_count, decreasing=TRUE)
    col_list = rownames(as.data.frame(col_count_sorted))[1:20]
  }
  if(missing(max_pro)){
    max_pro = ceiling(10*max(data))/10
  }
  data1 = data[1:46,col_list]
  data2 = data[47:92,col_list]
  if(clean==TRUE){
    colnames(data1) = as.character(sapply(colnames(data1), clean_met_name))
    colnames(data2) = as.character(sapply(colnames(data2), clean_met_name))
  }
  pm1 = pheatmap(data1, cluster_rows=FALSE, cluster_cols=FALSE,
                 border=FALSE,breaks=seq(0,max_pro,length.out=100),
                 cellwidth = 15, cellheight =9, fontsize=9, legend=FALSE)
  pm2 = pheatmap(data2, cluster_rows=FALSE, cluster_cols=FALSE,
                 border=FALSE,breaks=seq(0,max_pro,length.out=100),
                 cellwidth = 15, cellheight =9, fontsize=9)
  grid.arrange(grobs=list(pm1[[4]],pm2[[4]]),ncol=2)
}

clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

# data
metab = readRDS('../../used/metab_1776.rds')
metab = as.matrix(metab)
prot = readRDS('../../used/proteins_denoised_new.rds')
cov = readRDS('../../used/covars_587.rds')

metab_ori_neg = read_csv('../../Original_data/Lung cancer feature table RP NEG_231019.csv')
metab_ori_pos = read_csv('../../Original_data/Lung cancer feature table RP POS_231019.csv')

sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever.rds'))
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))

sel_pro_ever_subsample = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever_subsample.rds'))

dim(sel_pro_ever_subsample)

# compare sel_pro_ever with sel_pro_ever_subsample


plot_2_heatmap(sel_pro_ever,clean=TRUE)
plot_2_heatmap(sel_pro_ever_subsample,clean=TRUE)


sel_pro_ever_count_prot = sort(apply(sel_pro_ever, 1, function(x){sum(x!=0)}), decreasing=TRUE)
rownames(as.data.frame(sel_pro_ever_count_prot))

# Sort by proteins
sel_pro_ever_sorted = sel_pro_ever[rownames(as.data.frame(sel_pro_ever_count_prot)),]
sel_pro_ever_subsample_sorted = sel_pro_ever_subsample[rownames(as.data.frame(sel_pro_ever_count_prot)),]

plot_2_heatmap(sel_pro_ever_sorted,clean=TRUE)
plot_2_heatmap(sel_pro_ever_subsample_sorted,clean=TRUE)

# Sort by most selected metablites in whole dataset
met_list_ever_20 = apply(sel_pro_ever,2,function(x){sum(x!=0)})
met_list_ever_20 = sort(met_list_ever_20,decreasing=TRUE)[1:20]
met_list_ever_20 = rownames(as.data.frame(met_list_ever_20))


plot_2_heatmap(sel_pro_ever_sorted, col_list = met_list_ever_20, clean=TRUE)
plot_2_heatmap(sel_pro_ever_subsample_sorted, col_list = met_list_ever_20, clean=TRUE)
