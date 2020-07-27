path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))

suppressPackageStartupMessages(library(mixOmics))

StabilityPlot=function(X, Y, NIter=100, plot=FALSE){
  MyStab=NULL
  
  pb=txtProgressBar(style=3)
  for (NVar in seq((ncol(Y)-1), 1, by=-1)){
    setTxtProgressBar(pb, 1-(NVar-1)/(ncol(X)))
    TmpStab=NULL
    for (k in 1:NIter){
      s=sample(seq(1, nrow(X)), size=0.8*nrow(X), replace = FALSE) # subsampling procedure
      X_sub=X[s,]
      rownames(X_sub)=s
      Y_sub=Y[s,]
      
      TmpsPLS=spls(X_sub, Y_sub, keepY=NVar, ncomp=1, mode='regression')
      TmpStab=rbind(TmpStab, TmpsPLS$loadings$Y[,1])
    }
    MyStab=rbind(MyStab, (apply(TmpStab, 2, FUN=function(x){sum(x!=0)})/NIter))
  }
  
  MyStab=rbind(rep(1, ncol(Y)), MyStab)
  rownames(MyStab)=seq(ncol(Y), 1, by=-1)
  
  if (plot){
    pheatmap(MyStab, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
  }
  
  return(MyStab)
}

set.seed(1)
prot_stability = StabilityPlot(X = metab, Y = prot, NIter = 100)


saveRDS(prot_stability,
        file = "/rdsgpfs/general/project/hda_tds/live/GROUP_2/results/Lucas/sPLS/prot_stability.rds")
