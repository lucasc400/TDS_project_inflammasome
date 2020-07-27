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
met_sel_count_df=as.data.frame(met_sel_count)
pro_sel_count = apply(lasso_beta, 1, function(x){sum(x!=0)} )

# Cleaning metabolites name
clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

# LASSO -------------------------------------------------------------------
# Summarise LASSO results
sum(met_sel_count!=0)
top20 = as.data.frame(met_sel_count[1:20])
#rownames(top20) = as.character(sapply(rownames(top20), clean_met_name))
top20_names = rownames(top20)

# Table for stratified LASSO results
beta_cancer = readRDS('beta_cancer.rds')
beta_smoking = readRDS('beta_smoking.rds')

beta_case = beta_cancer$case
beta_control = beta_cancer$control
beta_ever = beta_smoking$ever
beta_never = beta_smoking$never

beta_case_count = apply(beta_case,2,function(x){sum(x!=0)})
beta_control_count = apply(beta_control,2,function(x){sum(x!=0)})
beta_ever_count = apply(beta_ever,2,function(x){sum(x!=0)})
beta_never_count = apply(beta_never,2,function(x){sum(x!=0)})

beta_case_count_pro = apply(beta_case,1,function(x){sum(x!=0)})
beta_control_count_pro = apply(beta_control,1,function(x){sum(x!=0)})
beta_ever_count_pro = apply(beta_ever,1,function(x){sum(x!=0)})
beta_never_count_pro = apply(beta_never,1,function(x){sum(x!=0)})


create_beta_matrix <- function(){
  M = matrix(0,dim(prot)[2],dim(metab)[2])
  colnames(M) = colnames(metab)
  rownames(M) = colnames(prot)
  return(M)
}

beta_case_logical = (beta_case!=0)
beta_control_logical = (beta_control!=0)

beta_cancer_diff = (beta_case_logical != beta_control_logical)
beta_cancer_diff = apply(beta_cancer_diff,2,as.numeric)

beta_cancer_int = (beta_case_logical == 1 & beta_control_logical == 1)
beta_cancer_int = apply(beta_cancer_int,2,as.numeric)

beta_ever_logical = (beta_ever!=0)
beta_never_logical = (beta_never!=0)

beta_smoking_diff = (beta_ever_logical != beta_never_logical)
beta_smoking_diff = apply(beta_smoking_diff,2,as.numeric)

beta_smoking_int = (beta_ever_logical == 1 & beta_never_logical == 1)
beta_smoking_int = apply(beta_smoking_int,2,as.numeric)

beta_cancer_diff_count = apply(beta_cancer_diff,2,function(x){sum(x!=0)})
beta_cancer_int_count = apply(beta_cancer_int,2,function(x){sum(x!=0)})
beta_smoking_diff_count = apply(beta_smoking_diff,2,function(x){sum(x!=0)})
beta_smoking_int_count = apply(beta_smoking_int,2,function(x){sum(x!=0)})

cancer_table <- data.frame('All'=met_sel_count[top20_names],
                           'Case'=beta_case_count[top20_names],
                           'Control'=beta_control_count[top20_names],
                           'Differential'=beta_cancer_diff_count[top20_names],
                           'Intersection'=beta_cancer_int_count[top20_names])
rownames(cancer_table) = as.character(sapply(rownames(cancer_table),clean_met_name))
cancer_table

smoking_table <- data.frame('All'=met_sel_count[top20_names],
                            'Ever smoker'=beta_ever_count[top20_names],
                            'Never smoker'=beta_never_count[top20_names],
                            'Differential'=beta_smoking_diff_count[top20_names],
                            'Intersection'=beta_smoking_int_count[top20_names])
rownames(smoking_table) = as.character(sapply(rownames(smoking_table),clean_met_name))
smoking_table


# Stability ---------------------------------------------------------------
sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_control = readRDS(paste0(path,'/Stability/Cancer/sel_pro_control.rds'))
sel_pro_case = readRDS(paste0(path,'/Stability/Cancer/sel_pro_case.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever.rds'))
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))

colnames(sel_pro_control) = as.character(sapply(colnames(sel_pro_control), clean_met_name))
colnames(sel_pro_case) = as.character(sapply(colnames(sel_pro_case), clean_met_name))
colnames(sel_pro_never) = as.character(sapply(colnames(sel_pro_never), clean_met_name))
colnames(sel_pro_ever) = as.character(sapply(colnames(sel_pro_ever), clean_met_name))

# Counts
met_sel_pro_count = apply(sel_pro, 2, function(x){sum(x!=0)})
pro_sel_pro_count = apply(sel_pro, 1, function(x){sum(x!=0)})

# Visualisation & heatmap
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

# Stratified heatmaps
met_sel_pro_sum = apply(sel_pro, 2, sum)
length(met_sel_pro_sum)
top20_pro_sum = sort(met_sel_pro_sum, decreasing=TRUE)[1:20]
top20_pro_sum_names = rownames(as.data.frame(top20_pro_sum))

plot_2_heatmap(sel_pro_control)
plot_2_heatmap(sel_pro_control,col_list = top20_pro_sum_names)

plot_2_heatmap(sel_pro_case)
plot_2_heatmap(sel_pro_case,col_list = top20_pro_sum_names)

plot_2_heatmap(sel_pro_ever)
plot_2_heatmap(sel_pro_ever,col_list = top20_pro_sum_names)

plot_2_heatmap(sel_pro_never)
plot_2_heatmap(sel_pro_never,col_list = top20_pro_sum_names)


# Analysis of Stability ---------------------------------------------------
# top 20 selected metabolites in stability
top20_pro = as.data.frame(sort(met_sel_pro_count,decreasing=TRUE)[1:20])
top20_pro
top20_pro_names = rownames(top20_pro)

# compare top 20 in LASSO and stability
top20_names_clean = as.character(sapply(top20_names, clean_met_name))
sum(top20_names_clean %in% top20_pro_names) # 11


lasso_stab_df = as.data.frame(sort(met_sel_count,decreasing=TRUE))
colnames(lasso_stab_df) = '#LASSO_selection'
rownames(lasso_stab_df) = as.character(sapply(rownames(lasso_stab_df), clean_met_name))

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

#lasso_stab_df[1:30,]
lasso_stab_df['232.8163@0.6834795',]

# Can I sort it by stability?
lasso_stab_df_new <- lasso_stab_df %>%
  rownames_to_column('metab_name') %>%
  arrange(desc(Stability.0.5)) %>%
  column_to_rownames('metab_name')

lasso_stab_df_new[1:30,]

# Proteins that had most and no selection
sort(pro_sel_pro_count,decreasing=TRUE)
pro_sel_pro_count[pro_sel_pro_count==0]

# Analysis of stratified LASSO --------------------------------------------
# Compare the most selected metabolites

# Cancer
beta_control_top20 = sort(beta_control_count,decreasing=TRUE)[1:20]
beta_case_top20 = sort(beta_case_count,decreasing=TRUE)[1:20]
beta_ever_top20 = sort(beta_ever_count,decreasing=TRUE)[1:20]
beta_never_top20 = sort(beta_never_count,decreasing=TRUE)[1:20]

beta_control_top20 %in% beta_case_top20

prot_cancer_df = as.data.frame(sort(beta_case_count_pro,decreasing=TRUE))
prot_cancer_df$control = as.numeric(as.data.frame(sort(beta_control_count_pro,decreasing=TRUE))[rownames(prot_cancer_df),])
colnames(prot_cancer_df)[1] = 'case'

prot_cancer_df_nozero = prot_cancer_df[apply(prot_cancer_df, 1, function(x){any(x!=0)}),]
ggplot(data=prot_cancer_df_nozero, aes(x=control,y=case,label=name))+geom_point()+
  geom_text(aes(label=name),hjust=0,vjust=0)+
  geom_abline(slope=1,intercept=0)

prot_cancer_df_nozero

# Compare the proteins with most selections, and those with no selections


#met_sel_pro_count = apply(sel_pro, 2, function(x){sum(x!=0)})

# how similar are LASSO and stability
met_sel_count_clean = as.data.frame(met_sel_count)
rownames(met_sel_count_clean) = sapply(rownames(met_sel_count_clean), clean_met_name)
mel_sel_count_clean = as.vector(met_sel_count_clean)
sum(rownames(as.data.frame(sort(met_sel_pro_count,decreasing=TRUE)[1:50])) %in% rownames(met_sel_count_clean)[1:50])



# Analysis of stratified stability ----------------------------------------
stratified_list = list('control'=sel_pro_control, 'case'=sel_pro_case, 
                       'ever'=sel_pro_ever, 'never'=sel_pro_never)
sapply(stratified_list, function(x){sum(x!=0)})
sapply(stratified_list, function(x){sum(x>=0.5)})
sapply(stratified_list, function(x){sum(x>=0.6)})
sapply(stratified_list, function(x){sum(x>=0.7)})
sapply(stratified_list, function(x){sum(x>=0.8)})
sapply(stratified_list, function(x){sum(x>=0.9)})