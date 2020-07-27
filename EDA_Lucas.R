rm(list=ls())

# Working directories
path <- "/rdsgpfs/general/project/hda_tds/live/GROUP_2"
setwd(path)

# Libraries
library(pacman)
p_load(ggplot2, tidyverse, dplyr, GGally, rstudioapi, stringi,
       naniar, openxlsx)

# Functions
transpose_rowname <- function(X){
  X = as.data.frame(t(as.matrix(X)))
  colnames(X) = X[1,]
  X = X[2:nrow(X),]
  return(X)
}

# Quantile score ----------------------------------------------------------
protein = readRDS(paste0(path,'/used/proteins_denoised_new.rds'))
dim(protein)

# quantile(unlist(protein_clean),probs=0.75,na.rm=TRUE) # overall quantile
compute_qs <- function(X){
  stopifnot(any(class(X)=='data.frame'))
  p_3rdq = apply(X, 2, function(x) quantile(x,probs=0.75,na.rm=TRUE))
  score_matrix = ifelse(X>=p_3rdq,1,0)
  scores = apply(score_matrix,1,function(x) sum(x,na.rm=TRUE))
  return(scores)
}

qt_score = compute_qs(protein)
# covariates$qt_score = qt_score

# Density plot of qs scores
# ggplot(aes(x=qt_score),data=covariates)+geom_density()


# Output
write_rds(qt_score,paste0(path,'/used/proteins_qt_score_new.rds'))



# To be modified ----------------------------------------------------------



covariates <- read_csv(paste0(getwd(),"/proteins_covariates.csv"))
covariates <- covariates[,-1]
#proteins<-readRDS(paste0(getwd(),"/Proteins/Proteins.rds"))
#proteins_cov <- readRDS(paste0(getwd(),"/Proteins/Proteins_technical_covariates.rds"))
#metab_pos<-read_csv(paste0(getwd(),"/Metabolites/Lung cancer feature table RP POS_231019.csv"))
#metab_neg<-read_csv(paste0(getwd(),"/Metabolites/Lung cancer feature table RP NEG_231019.csv"))

summary(covariates)
head(covariates)



# NLR ---------------------------------------------------------------------
compute_NLR <- function(X){
  stopifnot(any(class(X)=='data.frame'))
  attach(X)
  lympho = B+CD4T+CD8T+NK
  NLR = Neutrophils/lympho
  detach(X)
  return(NLR)
}

covariates$NLR = compute_NLR(covariates)
# Note that NA: 72 76 137 230 236 274 313 375 433 533 575 638


# Plate -------------------------------------------------------------------

metabolites = read.csv('metabolites_10.csv')
protein_cov = readRDS(paste0(path,'/Original_data/Proteins_technical_covariates.rds'))

met_neg_cov = read.xlsx(paste0(path,'/Original_data/LC-MS_worklist_LunCan_RPneg.xlsx'))

dim(met_neg_cov)
length(unique(met_neg_cov$Position))
length(unique(protein_cov$Plate.ID))

# Metabolites -------------------------------------------------------------

# Transposing metabolites
metab_pos_t = transpose_rowname(metab_pos)
metab_neg_t = transpose_rowname(metab_neg)

library(stringi)
# Covariates - 648, Metabolites - 666, Proteins - 702.
# Repeats in protein
str_detect(rownames(proteins),'._2') # sum = 56, 56 people had two samples

cov_rn = rownames(covariates)
length(cov_rn)
p_rn = rownames(proteins)
p_rn_no2 = p_rn[!str_detect(p_rn,'._2')]
length(p_rn_no2)

stopifnot(sort(cov_rn)==sort(p_rn_no2)) # matched

# Repeats in metabolites
sum(str_detect(rownames(metab_neg_t),'.Blank.')) # 16 - 16 blanks
sum(str_detect(rownames(metab_pos_t),'.Blank.'))

m_neg_rn = rownames(metab_neg_t)
m_neg_rn_nb = m_neg_rn[!str_detect(m_neg_rn,'.Blank.')]
length(m_neg_rn_nb) # 650

# Id match
id_match_table = table(covariates$id_match)
id_match_table[id_match_table == 1]

# NA
# prop_miss(metab_neg_t)
# metab_neg_t_df = as.data.frame(metab_neg_t)
# metab_neg_na <- metab_neg_t_df %>%
#   miss_var_summary()
# 
# ggplot(aes(x=pct_miss),data=metab_neg_na) + geom_density()
# summary(metab_neg_na)