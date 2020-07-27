rm(list=ls())
#path=dirname(rstudioapi::getActiveDocumentContext()$path)
#path=str_remove(path,'/Scripts/.*')
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2'
setwd(path)
path

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(omics))

cov_met = read_rds('metabolites2.rds')

colnames(cov_met)[1:10]
sum(is.na(cov_met))
met_index = c(seq(7,ncol(cov_met)))
# 2103767 2106943``

met = cov_met[,met_index]
cov = cov_met[,1:6]
colnames(cov)

log_met = log(met)

f_org = mlmer(met~(1|Plate),data=cov,save.residuals=TRUE,save.ranks=FALSE)

f_log = mlmer(log_met~(1|Plate)+(1|Position),data=cov,save.residuals=TRUE,save.ranks=FALSE)

#dim(f0$residuals)
denoised_met = f_org$residuals
denoised_met_df = as.data.frame(denoised_met)

denoised_log_met = f_log$residuals
denoised_log_met_df = as.data.frame(denoised_log_met)


means_org = colMeans(met)
means_denoised = colMeans(denoised_med_df)
summary(means_denoised)

par(mfrow=c(1,1))
plot(log(means_org), log(means_denoised))
