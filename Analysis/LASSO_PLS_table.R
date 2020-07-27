rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas'
setwd(path)
# libraries
library(pacman)
p_load(pheatmap,grid,gridExtra,stringr,ggplot2, dplyr, tibble)

# data
metab = readRDS('../../used/metab_1776.rds')
metab = as.matrix(metab)
prot = readRDS('../../used/proteins_denoised_new.rds')
cov = readRDS('../../used/covars_587.rds')

metab_ori_neg = read.csv('../../Original_data/Lung cancer feature table RP NEG_231019.csv')
metab_ori_pos = read.csv('../../Original_data/Lung cancer feature table RP POS_231019.csv')

sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever_subsample.rds')) # subsampled results
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))

load('/rdsgpfs/general/project/hda_tds/live/GROUP_2/TDS_project/Scripts/Vic/HPC_jobs/PLS_new/more_pls/Creating_heatmaps/results_whole_pop/Selection_whole_pop.RData')
load('/rdsgpfs/general/project/hda_tds/live/GROUP_2/TDS_project/Scripts/Vic/HPC_jobs/PLS_new/more_pls/Creating_heatmaps/results_whole_pop/selection_stratified.RData')

# function
clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

# cleaning
rownames(selection_whole_pop) = sapply(rownames(selection_whole_pop), clean_met_name)
rownames(selection_stratified) = sapply(rownames(selection_stratified), clean_met_name)

# stability_ranking
sel_pro_count.5 = apply(sel_pro,2,function(x){sum(x>=0.5)})
sel_pro_count.5 = as.data.frame(sel_pro_count.5)
colnames(sel_pro_count.5) = 'count'


sel_pro_count.0 = apply(sel_pro,2,function(x){sum(x>0)})
sel_pro_count.0 = as.data.frame(sel_pro_count.0)
colnames(sel_pro_count.0) = 'count'

selection_whole_pop$type = rownames(selection_whole_pop) %in% metab_ori_neg$Compound
selection_whole_pop$type = ifelse(selection_whole_pop$type, 'Neg', 'Pos')

selection_whole_pop$Stability.0.count = sel_pro_count.0[rownames(selection_whole_pop),'count']
selection_whole_pop$Stability.5.count = sel_pro_count.5[rownames(selection_whole_pop),'count']

selection_whole_pop$Stability.0.rank = rank(-selection_whole_pop$Stability.0.count,ties.method='min')
selection_whole_pop$Stability.5.rank = rank(-selection_whole_pop$Stability.5.count,ties.method='min')


selection_whole_pop = selection_whole_pop[,c(1,2,4,6,5,7,3)]
selection_whole_pop[1:20,]

selection_whole_pop$Stability.0 = paste0(selection_whole_pop$Stability.0.count, 
                                         ' (', selection_whole_pop$Stability.0.rank, ')')

selection_whole_pop$Stability.5 = paste0(selection_whole_pop$Stability.5.count, 
                                         ' (', selection_whole_pop$Stability.5.rank, ')')

table1 = selection_whole_pop[,c(1,2,9)]
colnames(table1)[3] = c('LASSO selection probability>=0.5 (rank)')

table1[1:10,]
#save(table_concise,file='/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/LASSO_pls_table.RData')

# stratified
#selection_stratified
sel_pro_ever_count.5 = apply(sel_pro_ever,2,function(x){sum(x>=0.5)})
sel_pro_ever_count.5 = as.data.frame(sel_pro_ever_count.5)
colnames(sel_pro_ever_count.5) = 'count'
rownames(sel_pro_ever_count.5) = sapply(rownames(sel_pro_ever_count.5),clean_met_name)

sel_pro_never_count.5 = apply(sel_pro_never,2,function(x){sum(x>=0.5)})
sel_pro_never_count.5 = as.data.frame(sel_pro_never_count.5)
colnames(sel_pro_never_count.5) = 'count'
rownames(sel_pro_never_count.5) = sapply(rownames(sel_pro_never_count.5),clean_met_name)

sel_pro_ever_count.0 = apply(sel_pro_ever,2,function(x){sum(x>0)})
sel_pro_ever_count.0 = as.data.frame(sel_pro_ever_count.0)
colnames(sel_pro_ever_count.0) = 'count'
rownames(sel_pro_ever_count.0) = sapply(rownames(sel_pro_ever_count.0),clean_met_name)

sel_pro_never_count.0 = apply(sel_pro_never,2,function(x){sum(x>0)})
sel_pro_never_count.0 = as.data.frame(sel_pro_never_count.0)
colnames(sel_pro_never_count.0) = 'count'
rownames(sel_pro_never_count.0) = sapply(rownames(sel_pro_never_count.0),clean_met_name)


selection_stratified$Stability.5.ever.count = sel_pro_ever_count.5[rownames(selection_stratified),'count']
selection_stratified$Stability.5.ever.rank = rank(-selection_stratified$Stability.5.ever.count,ties.method='min')
selection_stratified$Stability.5.ever = paste0(selection_stratified$Stability.5.ever.count,
                                               ' (', selection_stratified$Stability.5.ever.rank, ')')

selection_stratified$Stability.5.never.count = sel_pro_never_count.5[rownames(selection_stratified),'count']
selection_stratified$Stability.5.never.rank = rank(-selection_stratified$Stability.5.never.count,ties.method='min')
selection_stratified$Stability.5.never = paste0(selection_stratified$Stability.5.never.count,
                                                ' (', selection_stratified$Stability.5.never.rank, ')')

selection_stratified$Stability.0.ever.count = sel_pro_ever_count.0[rownames(selection_stratified),'count']
selection_stratified$Stability.0.ever.rank = rank(-selection_stratified$Stability.0.ever.count,ties.method='min')
selection_stratified$Stability.0.ever = paste0(selection_stratified$Stability.0.ever.count,
                                               ' (', selection_stratified$Stability.0.ever.rank, ')')

selection_stratified$Stability.0.never.count = sel_pro_never_count.0[rownames(selection_stratified),'count']
selection_stratified$Stability.0.never.rank = rank(-selection_stratified$Stability.0.never.count,ties.method='min')
selection_stratified$Stability.0.never = paste0(selection_stratified$Stability.0.never.count,
                                                ' (', selection_stratified$Stability.0.never.rank, ')')

head(selection_stratified)

table2 = selection_stratified[,c(1:6,9,12,15,18)]
colnames(table2)[7:8] = c('Smoker LASSO selection probability>=0.5 (rank)','Non-smoker LASSO selection probability>=0.5 (rank)')
table2 = table2[,c(3,5,8,7)]

colnames(table2)
table2 <- table2 %>%
  rownames_to_column('metabolite') %>%
  arrange(comp1_NS) %>%
  column_to_rownames('metabolite')

table2[1:10,]


table1_sorted_lasso <- selection_whole_pop %>%
  rownames_to_column('metabolite') %>%
  select(Selection_C1, Stability.5.count, metabolite) %>%
  arrange(desc(Stability.5.count)) %>%
  column_to_rownames('metabolite')

head(table1_sorted_lasso, n=30)
#write.csv(table1, 'table1.csv')
#write.csv(table2, 'table2.csv')

table1_20 = table1[1:20,]
colnames(table1_20)[1] = 'sPLS selection'
table1_20=table1_20[,c(1,3)]
colnames(table1_20)[2] = 'LASSO selection proportion>=0.5 (rank)'


saveRDS(table1_20, '../LASSO_PLS_table.rds')


ggplot(data=selection_whole_pop, aes(x=Selection_C1,y=Stability.5.count))+
  geom_point()+geom_smooth()
