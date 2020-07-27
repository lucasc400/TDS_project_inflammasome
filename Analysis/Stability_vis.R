rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas'
setwd(path)

library(pacman)
p_load(pheatmap,grid,gridExtra,stringr,ggplot2, dplyr, tibble)

metab = readRDS('../../used/metab_1776.rds')
metab = as.matrix(metab)
prot = readRDS('../../used/proteins_denoised_new.rds')
cov = readRDS('../../used/covars_587.rds')

sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever.rds'))
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))
sel_pro_ever_subsample = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever_subsample.rds'))


# Counts
met_sel_pro_count = apply(sel_pro, 2, function(x){sum(x!=0)})
pro_sel_pro_count = apply(sel_pro, 1, function(x){sum(x!=0)})

# Cleaning metabolites names
clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

plot_2_heatmap <- function(data, clean=FALSE, top=20,col_list, max_pro){
  # Packages required: pheatmap, grid, gridExtra, stringr
  if(missing(col_list)){
    col_count = apply(data,2,sum)
    col_count_sorted = sort(col_count, decreasing=TRUE)
    col_list = rownames(as.data.frame(col_count_sorted))[1:top]
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


# Clean metabolite names
colnames(sel_pro_never) = as.character(sapply(colnames(sel_pro_never), clean_met_name))
colnames(sel_pro_ever) = as.character(sapply(colnames(sel_pro_ever), clean_met_name))
colnames(sel_pro_ever_subsample) = as.character(sapply(colnames(sel_pro_ever_subsample), clean_met_name))
colnames(sel_pro_qt) = as.character(sapply(colnames(sel_pro_qt), clean_met_name))


# Sorting proteins (rows)
pro_sel_pro_count_names = sort(pro_sel_pro_count,decreasing=TRUE)
pro_sel_pro_count_names = rownames(as.data.frame(pro_sel_pro_count_names))
sel_pro_sorted = sel_pro[pro_sel_pro_count_names,]
sel_pro_ever_sorted = sel_pro_ever[pro_sel_pro_count_names,]
sel_pro_ever_subsample_sorted = sel_pro_ever_subsample[pro_sel_pro_count_names,]
sel_pro_never_sorted = sel_pro_never[pro_sel_pro_count_names,]



# Plotting heatmaps
met_sel_pro_sum = apply(sel_pro, 2, sum)
length(met_sel_pro_sum)
top20_pro_sum = sort(met_sel_pro_sum, decreasing=TRUE)[1:20]
top20_pro_sum_names = rownames(as.data.frame(top20_pro_sum))

# Sort by most selected metablites in whole dataset
top_30_met = rownames(as.data.frame(sort(met_sel_pro_count,decreasing=TRUE)))[1:30]


plot_2_heatmap(sel_pro_sorted, col_list=top_30_met)
plot_2_heatmap(sel_pro_ever_sorted, col_list=top_30_met)
plot_2_heatmap(sel_pro_ever_subsample_sorted, col_list=top_30_met)
plot_2_heatmap(sel_pro_never_sorted, col_list=top_30_met)
