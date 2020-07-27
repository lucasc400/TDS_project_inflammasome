rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2'
setwd(path)

# Libraries
library(pacman)
p_load(ggplot2, tidyverse, dplyr, GGally, rstudioapi, stringi, naniar)

# Load data
metab_pos<-read_csv(paste0(getwd(),"/Old_data/meta_pos_t_clean.csv"))
metab_neg<-read_csv(paste0(getwd(),"/Old_data/meta_neg_t_clean.csv"))
metabolites<-read_csv(paste0(getwd(),"/Old_data/metabolites.csv"))

metab_pos_org<-read_csv(paste0(getwd(),
                               "/Original_data/Lung cancer feature table RP POS_231019.csv"))

metab_neg_org<-read_csv(paste0(getwd(),
                               "/Original_data/Lung cancer feature table RP NEG_231019.csv"))


# Keep metabolites with <=10% missing values ------------------------------
keep_10 <- function(df){
  df_10_name <- df %>%
    miss_var_summary() %>%
    filter(pct_miss<=10) %>%
    select(variable)
  df_10_name = as.vector(unlist(df_10_name))
  # Put sample IDs as first column
  ID_index = which(df_10_name=='X1')
  df_10_name = df_10_name[c(ID_index, seq(1,ID_index-1), seq(ID_index+1,length(df_10_name)))]
  df_10 = df[,df_10_name]
  return(df_10)
}

metab_pos_10 <- keep_10(metab_pos)
metab_neg_10 <- keep_10(metab_neg)
metab <- keep_10(metabolites)

View(metab_pos_10)
View(metab_neg_10)
View(metab)

# output
write_csv(metab_pos_10,'meta_pos_10.csv')
write_csv(metab_neg_10,'meta_neg_10.csv')
write_csv(metab,'metabolites_10.csv')




# Separating RT and m/z (not used yet) ------------------------------------

# remove the first repeat, and ':' in the second

colnames(metab_neg)[which(str_detect(colnames(metab_neg),':'))][2] <- 
           str_extract(colnames(metab_neg)[which(str_detect(colnames(metab_neg),':'))][2],
            '^([^:])+')

metab_neg[,-which(str_detect(colnames(metab_neg),':'))]
dim(metab_neg)


neg_m = colnames(metab_neg)
neg_m = neg_m[-1]
length(neg_m)

sum(str_detect(neg_m,'@'))

peaks = str_extract(neg_m,'^([^@])+') # Everything before @
rts = str_extract(neg_m, '(?<=@).*')
length(rts)
length(unique(rts))