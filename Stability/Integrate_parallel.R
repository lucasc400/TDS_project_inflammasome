rm(list=ls())
path = "/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas/Stability/"
setwd(path)
library(ggplot2)
suppressPackageStartupMessages(library(pheatmap))
library(stringr)

sel_pro.1 = readRDS('All/sel_pro_1_new.rds')
sel_pro.2 = readRDS('All/sel_pro_2_new.rds')
sel_pro.3 = readRDS('All/sel_pro_3_new.rds')
sel_pro.4 = readRDS('All/sel_pro_4_new.rds')

prot = readRDS('../../../used/proteins_denoised_new.rds')

rownames(sel_pro.1) = colnames(prot)[1:23]
rownames(sel_pro.2) = colnames(prot)[24:46]
rownames(sel_pro.3) = colnames(prot)[47:69]
rownames(sel_pro.4) = colnames(prot)[70:92]

sel_pro = rbind(sel_pro.1,sel_pro.2,sel_pro.3,sel_pro.4)

clean_met_name <- function(x){
  # Converts 'X1111.111.1.1111' into '1111.111@1.1111'
  # Package required: stringr
  x = str_replace_all(x,'X','')
  x_vec = strsplit(x,'\\.')
  x_vec = x_vec[[1]]
  return(paste0(x_vec[1],'.',x_vec[2],'@',x_vec[3],'.',x_vec[4]))
}

colnames(sel_pro) = as.character(sapply(colnames(sel_pro),clean_met_name))

rowsum1 = as.data.frame(apply(sel_pro,1,sum))
colnames(rowsum1)='RS'
ggplot(aes(x=RS),data=rowsum1)+geom_density()

plot(density(sel_pro))
hist(sel_pro)

saveRDS(sel_pro,'new_stability.rds')




# Cancer ------------------------------------------------------------------


sel_pro_case.1 = readRDS(paste0(path,'Cancer/sel_pro_case_1_new.rds'))
sel_pro_case.2 = readRDS(paste0(path,'Cancer/sel_pro_case_2_new.rds'))
sel_pro_case.3 = readRDS(paste0(path,'Cancer/sel_pro_case_3_new.rds'))
sel_pro_case.4 = readRDS(paste0(path,'Cancer/sel_pro_case_4_new.rds'))

sel_pro_case = rbind(sel_pro_case.1, sel_pro_case.2, sel_pro_case.3, sel_pro_case.4)
colnames(sel_pro_case) = as.character(sapply(colnames(sel_pro_case),clean_met_name))

sel_pro_control.1 = readRDS(paste0(path,'Cancer/sel_pro_control_1_new.rds'))
sel_pro_control.2 = readRDS(paste0(path,'Cancer/sel_pro_control_2_new.rds'))
sel_pro_control.3 = readRDS(paste0(path,'Cancer/sel_pro_control_3_new.rds'))
sel_pro_control.4 = readRDS(paste0(path,'Cancer/sel_pro_control_4_new.rds'))

sel_pro_control = rbind(sel_pro_control.1, sel_pro_control.2, sel_pro_control.3, sel_pro_control.4)
colnames(sel_pro_control) = as.character(sapply(colnames(sel_pro_control),clean_met_name))

sum(sel_pro_control!=0)
sum(sel_pro_case!=0)

#pheatmap(t(sel_pro_control.1), cluster_rows=FALSE, cluster_cols=FALSE,
#         border=FALSE,breaks=seq(0,1,length.out=100))
saveRDS(sel_pro_case, 'new_stability_case.rds')
saveRDS(sel_pro_control, 'new_stability_control.rds')

# Smoking -----------------------------------------------------------------
sel_pro_never.1 = readRDS(paste0(path,'Smoking/sel_pro_never_1_new.rds'))
sel_pro_never.2 = readRDS(paste0(path,'Smoking/sel_pro_never_2_new.rds'))
sel_pro_never.3 = readRDS(paste0(path,'Smoking/sel_pro_never_3_new.rds'))
sel_pro_never.4 = readRDS(paste0(path,'Smoking/sel_pro_never_4_new.rds'))

sel_pro_never = rbind(sel_pro_never.1, sel_pro_never.2, sel_pro_never.3, sel_pro_never.4)
colnames(sel_pro_never) = as.character(sapply(colnames(sel_pro_never),clean_met_name))

sel_pro_ever.1 = readRDS(paste0(path,'Smoking/sel_pro_ever_1_new.rds'))
sel_pro_ever.2 = readRDS(paste0(path,'Smoking/sel_pro_ever_2_new.rds'))
sel_pro_ever.3 = readRDS(paste0(path,'Smoking/sel_pro_ever_3_new.rds'))
sel_pro_ever.4 = readRDS(paste0(path,'Smoking/sel_pro_ever_4_new.rds'))

sel_pro_ever = rbind(sel_pro_ever.1, sel_pro_ever.2, sel_pro_ever.3, sel_pro_ever.4)
colnames(sel_pro_ever) = as.character(sapply(colnames(sel_pro_ever),clean_met_name))

sum(sel_pro_never!=0)
sum(sel_pro_ever!=0)

saveRDS(sel_pro_never, 'new_stability_never.rds')
saveRDS(sel_pro_ever, 'new_stability_ever.rds')
