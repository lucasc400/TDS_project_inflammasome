rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas'
setwd(path)

library(pacman)
p_load(pheatmap,grid,gridExtra,stringr,ggplot2, dplyr, tibble)

metab = readRDS('../../used/metab_1776.rds')
metab = as.matrix(metab)
prot = readRDS('../../used/proteins_denoised_new.rds')
cov = readRDS('../../used/covars_587.rds')
lasso_beta = readRDS('metab_beta_lasso.rds')
met_sel_count = readRDS('metab_sel_count.rds')

sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_control = readRDS(paste0(path,'/Stability/Cancer/sel_pro_control.rds'))
sel_pro_case = readRDS(paste0(path,'/Stability/Cancer/sel_pro_case.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever.rds'))
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))

# Counts
met_sel_pro_count = apply(sel_pro, 2, function(x){sum(x!=0)})
pro_sel_pro_count = apply(sel_pro, 1, function(x){sum(x!=0)})
# Cleaning metabolites name
clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

names(met_sel_count) = sapply(names(met_sel_count), clean_met_name)

# Analysis of Stability ---------------------------------------------------
# top 20 selected metabolites in stability
top20_pro = as.data.frame(sort(met_sel_pro_count,decreasing=TRUE)[1:20])
top20_pro
top20_pro_names = rownames(top20_pro)


lasso_stab_df = as.data.frame(sort(met_sel_count,decreasing=TRUE))
colnames(lasso_stab_df) = '#LASSO_selection'

lasso_stab_df$Stability.0 = as.vector(as.data.frame(met_sel_pro_count)[rownames(lasso_stab_df),])

met_sel_pro_count.0.5 = apply(sel_pro, 2, function(x){sum(x>=0.5)})
lasso_stab_df$Stability.0.5 = as.vector(as.data.frame(met_sel_pro_count.0.5)[rownames(lasso_stab_df),])

met_sel_pro_count.0.8 = apply(sel_pro, 2, function(x){sum(x>=0.8)})
lasso_stab_df$Stability.0.8 = as.vector(as.data.frame(met_sel_pro_count.0.8)[rownames(lasso_stab_df),])

met_sel_pro_count.0.9 = apply(sel_pro, 2, function(x){sum(x>=0.9)})
lasso_stab_df$Stability.0.9 = as.vector(as.data.frame(met_sel_pro_count.0.9)[rownames(lasso_stab_df),])

lasso_stab_df[1:30,]

# append selection probability with qt_score
sel_pro_qt_df = as.data.frame(sel_pro_qt)
colnames(sel_pro_qt_df) = sapply(colnames(sel_pro_qt_df), clean_met_name)
lasso_stab_df$Stability.qt = as.numeric(sel_pro_qt_df[,rownames(lasso_stab_df)])


# Sort by Stability.0.5
lasso_stab_df_new <- lasso_stab_df %>%
  rownames_to_column('metab_name') %>%
  arrange(desc(Stability.0.5)) %>%
  column_to_rownames('metab_name')

lasso_stab_df_new[1:30,]





# Analysis of stratified stability ----------------------------------------
stratified_list = list('control'=sel_pro_control, 'case'=sel_pro_case, 
                       'ever'=sel_pro_ever, 'never'=sel_pro_never)
sapply(stratified_list, function(x){sum(x!=0)})
sapply(stratified_list, function(x){sum(x>=0.5)})
sapply(stratified_list, function(x){sum(x>=0.8)})
sapply(stratified_list, function(x){sum(x>=0.9)})



# Proteins that had most and no selection
sort(pro_sel_pro_count,decreasing=TRUE)
pro_sel_pro_count[pro_sel_pro_count==0]

