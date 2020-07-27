# path
rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas'
setwd(path)

# package
library(pacman)
p_load(pheatmap,grid,gridExtra,stringr,ggplot2, dplyr, tibble, ggrepel)

# loading data
metab = readRDS('../../used/metab_1776.rds')
metab = as.matrix(metab)
prot = readRDS('../../used/proteins_denoised_new.rds')
cov = readRDS('../../used/covars_587.rds')

sel_pro = readRDS(paste0(path,'/Stability/new_stability.rds'))
sel_pro_control = readRDS(paste0(path,'/Stability/Cancer/sel_pro_control.rds'))
sel_pro_case = readRDS(paste0(path,'/Stability/Cancer/sel_pro_case.rds'))
sel_pro_ever = readRDS(paste0(path,'/Stability/Smoking/sel_pro_ever.rds'))
sel_pro_never = readRDS(paste0(path,'/Stability/Smoking/sel_pro_never.rds'))
sel_pro_qt = readRDS(paste0(path, '/Stability/qt_stability.rds'))

prot_stability = readRDS(paste0(path,'/sPLS/prot_stability.rds'))
prot_stability_sum = apply(prot_stability,2,sum)
prot_stability_sum = sort(prot_stability_sum, decreasing=TRUE)
prot_stability_sorted_list = rownames(as.data.frame(prot_stability_sum))

pheatmap(prot_stability[,prot_stability_sorted_list], cluster_rows=FALSE, cluster_cols=FALSE,
         fontsize=10,cellwidth=12,cellheight=9,filename='protein_heatmap.png')

dim(sel_pro)
sel_pro_count_prot = apply(sel_pro, 1, function(x){sum(x!=0)})
sort(sel_pro_count_prot)

prot_table = as.data.frame(prot_stability_sum)
prot_table$LASSO_stab_count.0 = as.data.frame(sel_pro_count_prot)[rownames(prot_table),]



for (i in (seq(1:92))){
  pname<-colnames(prot_stability)[i]
  j<-1
  while (j<92 && prot_stability[j,pname]>0.5){
    j<-j+1
  }
  prot_table[pname,'pls_selection'] = 92-(j-1)
  
}

prot_table[,2:3] %>%
  rownames_to_column('prot_name') %>%
  arrange(desc(pls_selection)) %>%
  column_to_rownames('prot_name')

cor(prot_table$pls_selection, prot_table$LASSO_stab_count.0)

prot_table$name = rownames(prot_table)
colnames(prot_table)[3] = 'sPLS_selection'

ggplot(data=prot_table, aes(x=sPLS_selection,y=LASSO_stab_count.0,label=name)) + geom_point(colour='red') +
  geom_smooth() + 
  geom_text_repel(aes(label=name),hjust=0,vjust=0,size=5) + 
  scale_x_continuous(name='Minimum number of proteins needed in sPLS to reach 0.5 selection proportion') +
  scale_y_continuous(name='Number of metabolites selected in LASSO stability') +
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))
ggsave('LASSO_pls_point.png',width=55,height=35,unit='cm')
