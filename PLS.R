rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)

metab = readRDS(paste0(path,'metab_denoised_definite.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score.rds'))

suppressPackageStartupMessages(library(mixOmics))


# Tuning sPLS
t0 = Sys.time()
spls_tuning = tune.spls(metab, prot, test.keepX=c(1:200), progressBar=TRUE, nrepeat=5)
t1 = Sys.time()
t1-t0

spls_tuning$choice.keepX # optimal keepX=28

# sPLS
spls_model = spls(X=metab, Y=prot, keepX=28, mode='regression')
spls_model$explained_variance

# The 28 metabolites selected and their loading
metab_selected = as.data.frame(spls_model$loadings$X[,1][spls_model$loadings$X[,1]!=0])
#Note: spls_model$loadings return 'two components' for both X and Y - I'm not what that means
colnames(metab_selected)='loading'
metab_selected
