rm(list=ls())
library(sgPLS)
library(pheatmap)
library(parallel)

path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score_new.rds'))



#################### Functions ####################

SetToCoordinate=function(matrix, rown, coln, value){
  matrix[rown+nrow(matrix)*(coln-1)]=value
  return(matrix)
}


PlotCalibsPLS=function(res, sparsity='X', ncomp=1, main=NULL){
  if (sparsity%in%c('X', 'Y')){
    plot(res$MSEP[res$NComp==ncomp],
         type="l", xaxt="n", xlab=paste0("Number of variables in ", sparsity), ylab="MSEP",
         main=ifelse(is.null(main), yes=paste0("Calibration of the number of variables in ", sparsity), no=main),
         sub=paste("Number of components =", ncomp), las=1)
    axis(side = 1, at = 1:sum(res$NComp==ncomp))
    points(which.min(res$MSEP[res$NComp==ncomp]),
           res$MSEP[res$NComp==ncomp][which.min(res$MSEP[res$NComp==ncomp])],
           pch=19, col="tomato")
  } else {
    Tmp=matrix(0, nrow=length(unique(res$NVarX)), ncol=length(unique(res$NVarY)))
    Tmp=SetToCoordinate(Tmp, res$NVarX, res$NVarY, res$MSEP)
    
    rownames(Tmp)=unique(res$NVarX)
    colnames(Tmp)=unique(res$NVarY)
    pheatmap(Tmp, cluster_rows = F, cluster_cols = F)
  }
}


CalibratesPLS<-function(X, Y, nComp=1, nrepeat=1000, sparsity='X', 
                        plot=TRUE, ask=FALSE, stepX=1, stepY=1,
                        MaxVarX=ncol(X), MaxVarY=ncol(Y)){
  # plot: logical, used only if nComp is set to 1.
  # If TRUE, the calibration plot is done
  
  # ask: logical, used only if nComp is set to 1 and plot to TRUE.
  # If TRUE, the user is asked for the number of variables to include in the model
  # and this number is returned in the output
  
  # To be used locally if nComp>1 
  
  if(!sparsity%in%c('X', 'Y', 'XY')){
    stop('sparsity should be one of "X", "Y" and "XY"')
  }
  
  SelectedX<-NULL
  SelectedY<-NULL
  
  MinVarX=ifelse(sparsity%in%c('X', 'XY'), yes=1, no=MaxVarX)
  MinVarY=ifelse(sparsity%in%c('Y', 'XY'), yes=1, no=MaxVarY)
  
  Summary <- data.frame(NComp=NULL, NVarX=NULL, NVarY=NULL, MSEP=NULL, Q2=NULL)
  
  for(comp in c(1:nComp)){
    print(comp)
    
    if (comp>1){
      PlotCalibsPLS(Summary, ncomp=comp-1)
      print(Summary)
      
      if (sparsity%in%c('X', 'XY')){
        previousX=0
        while(!previousX%in%seq(1,MaxVarX)){
          previousX=as.numeric(as.character(readline(paste0('Number of variables in X in component ', comp-1,
                                                            '? (give a value between 1 and ', MaxVarX, ') '))))
        }
        SelectedX=c(SelectedX, previousX)
      }
      
      if (sparsity%in%c('Y', 'XY')){
        previousY=0
        while(!previousY%in%seq(1,MaxVarY)){
          previousY=as.numeric(as.character(readline(paste0('Number of variables in Y in component ', comp-1,
                                                            '? (give a value between 1 and ', MaxVarY, ') '))))
        }
        SelectedY=c(SelectedY, previousY)
      }
    }
    
    no_cores <- detectCores()-1
    clust <- makeCluster(no_cores, type='FORK')
    
    v=expand.grid(seq(MinVarY, MaxVarY, by=stepY), seq(MinVarX, MaxVarX, by=stepX))
    Summary=t(parApply(clust, v, MARGIN=1, FUN=function(k){
      j=k[1];i=k[2]
      TmpKeepX <- c(SelectedX,i)
      TmpKeepY <- c(SelectedY,j)
      TmpsPLS <- spls(X,Y,keepX=TmpKeepX, keepY=TmpKeepY,ncomp=comp,mode='regression')
      TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=nrepeat,progressBar=FALSE)
      return(c(comp, i, j, apply(TmpPerf$MSEP, 2, mean)[comp], TmpPerf$Q2.total[comp,]))
    }))
    
    stopCluster(clust)
    Summary=as.data.frame(Summary)
    colnames(Summary)=c('NComp', 'NVarX', 'NVarY', 'MSEP', 'Q2')
  }
  
  if (plot|(nComp>1)){
    PlotCalibsPLS(Summary, sparsity = sparsity, ncomp=nComp)
    print(Summary)
    
    if (ask|(nComp>1)){
      if (sparsity%in%c('X', 'XY')){
        previousX=0
        while(!previousX%in%seq(1,MaxVarX)){
          previousX=as.numeric(as.character(readline(paste0('Number of variables in X in component ', comp,
                                                            '? (give a value between 1 and ', MaxVarX, ') '))))
        }
        SelectedX=c(SelectedX, previousX)
      }
      
      if (sparsity%in%c('Y', 'XY')){
        previousY=0
        while(!previousY%in%seq(1,MaxVarY)){
          previousY=as.numeric(as.character(readline(paste0('Number of variables in Y in component ', comp,
                                                            '? (give a value between 1 and ', MaxVarY, ') '))))
        }
        SelectedY=c(SelectedY, previousY)
      }
    }
  }
  
  if (nComp>1|(plot&ask)){
    if (sparsity=='X'){
      return(list(Summary=Summary, keepX=SelectedX))
    }
    if (sparsity=='Y'){
      return(list(Summary=Summary, keepY=SelectedY))
    }
    if (sparsity=='XY'){
      return(list(Summary=Summary, keepX=SelectedX, keepY=SelectedY))
    }
  } else {
    return(Summary)
  }
}



CalibratesPLSperComp<-function(X, Y, comp=1, SelectedX=NULL, SelectedY=NULL, nrepeat=1000, sparsity='X', plot=TRUE, ask=FALSE){
  # SelectedX: number of variables in X selected for previous components, used only if sparsity is X or XY
  # SelectedY: number of variables in X selected for previous components, used only if sparsity is Y or XY  
  # Both SelectedX and SelectedY are used only if comp>1
  
  if(!sparsity%in%c('X', 'Y', 'XY')){
    stop('sparsity should be one of "X", "Y" and "XY"')
  }
  
  MaxVarX <- ncol(X)
  MaxVarY <- ncol(Y)
  
  MinVarX=ifelse(sparsity%in%c('X', 'XY'), yes=1, no=MaxVarX)
  MinVarY=ifelse(sparsity%in%c('Y', 'XY'), yes=1, no=MaxVarY)
  
  Summary <- data.frame(NComp=NULL, NVarX=NULL, NVarY=NULL, MSEP=NULL, Q2=NULL)
  
  print(comp)
  
  no_cores <- detectCores()-1
  clust <- makeCluster(no_cores, type='FORK')
  
  v=expand.grid(MinVarY:MaxVarY, MinVarX:MaxVarX)
  Summary=t(parApply(clust, v, MARGIN=1, FUN=function(k){
    j=k[1];i=k[2]
    TmpKeepX <- c(SelectedX,i)
    TmpKeepY <- c(SelectedY,j)
    TmpsPLS <- spls(X,Y,keepX=TmpKeepX, keepY=TmpKeepY,ncomp=comp,mode='regression')
    TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=nrepeat,progressBar=FALSE)
    return(c(comp, i, j, apply(TmpPerf$MSEP, 2, mean)[comp], TmpPerf$Q2.total[comp,]))
  }))
  
  stopCluster(clust)
  Summary=as.data.frame(Summary)
  colnames(Summary)=c('NComp', 'NVarX', 'NVarY', 'MSEP', 'Q2')
  
  return(Summary)
}



# Toy codes ---------------------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene[,1:20]
Y <- liver.toxicity$clinic
##### Calibration

set.seed(1)
res=CalibratesPLS(X, Y, nComp = 1, nrepeat=1, sparsity='XY')
## This time the minimal MSEP values is assessed on a grid of 
## number of variables in X and Y
res[which.min(res$MSEP),]
## In this example, I chose 19 variables in 
## and 5 variables in Y

# Analysis ----------------------------------------------------------------

# One component, sparsity on X
t0 = Sys.time()
set.seed(1)
res_x=CalibratesPLS(metab, prot, nComp = 1, nrepeat=1, sparsity='X')
PlotCalibsPLS(res_x, sparsity='X', ncomp=1) # calibration plot
t1 = Sys.time()
print(t1-t0)

# One component, sparsity on Y
set.seed(1)
res_y=CalibratesPLS(metab, prot, nComp = 1, nrepeat=1, sparsity='Y')
PlotCalibsPLS(res_y, sparsity='Y', ncomp=1) # calibration plot

# One component, sparsity on X and Y
res_xy = CalibratesPLS(metab, prot, nComp=1, nrepeat=1, sparsity='XY', plot=FALSE, MaxVarX = 300, MaxVarY = 92)


saveRDS(res_xy,'../results/Lucas/sPLS/calibration_XY.rds')

dim(res_xy)
res_xy_MSEP = res_xy[, c('NVarX', 'NVarY', 'MSEP')]
res_xy_MSEP_wide = pivot_wider(res_xy_MSEP, names_from=NVarY, values_from=MSEP)
res_xy_MSEP_wide = res_xy_MSEP_wide[,-1]
dim(res_xy_MSEP_wide)
res_xy_MSEP_wide.matrix = t(as.matrix(res_xy_MSEP_wide))

pheatmap(res_xy_MSEP_wide.matrix,cluster_rows=FALSE, cluster_cols=FALSE)


t0 = Sys.time()
res_xy = CalibratesPLS(metab, prot, nComp=1, nrepeat=1, sparsity='XY', plot=FALSE, MaxVarX = 1776, MaxVarY = 92)
saveRDS(res_xy,'../results/Lucas/sPLS/calibration_XY_1776.rds')
t1 = Sys.time()
print(t1-t0)

res_xy_MSEP = res_xy[, c('NVarX', 'NVarY', 'MSEP')]
res_xy_MSEP_wide = pivot_wider(res_xy_MSEP, names_from=NVarY, values_from=MSEP)
res_xy_MSEP_wide = res_xy_MSEP_wide[,-1]
dim(res_xy_MSEP_wide)
res_xy_MSEP_wide.matrix = t(as.matrix(res_xy_MSEP_wide))
saveRDS(res_xy_MSEP_wide.matrix, '../results/Lucas/sPLS/calibration_XY_1776_MSEP.rds')
