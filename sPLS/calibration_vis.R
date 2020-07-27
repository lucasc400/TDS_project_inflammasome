rm(list=ls())

path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas/sPLS/'
setwd(path)

library(pheatmap)
library(ggplot2)
library(sgPLS)
library(mixOmics)
calibration_xy = readRDS('calibration_XY_1776.rds')
MSEP_xy = readRDS('calibration_XY_1776_MSEP.rds')
metab = readRDS(paste0(path, '../../../used/metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'../../../used/proteins_denoised_new.rds'))
cov = readRDS(paste0(path, '../../../used/covars_587.rds'))


pheatmap(MSEP_xy,cluster_rows=FALSE, cluster_cols=FALSE)

plot(x=1:1776, y=MSEP_xy[76,], type='l')

which(MSEP_xy==min(MSEP_xy),arr.ind=TRUE)


MysPLSXY <- spls(metab,prot,keepX=24, keepY=76, ncomp=1,mode='regression')

MysPLSXY$loadings$X[MysPLSXY$loadings$X!=0]
MysPLSXY$loadings$Y

plotLoadings(MysPLSXY)

colnames(calibration_xy)
calibration_xy[calibration_xy$NVarX==24 & calibration_xy$NVarY==76,]
