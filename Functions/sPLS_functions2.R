rm(list=ls())
library(sgPLS)
library(pheatmap)
library(parallel)

data(liver.toxicity)
X <- liver.toxicity$gene[,1:20]
Y <- liver.toxicity$clinic


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


CalibratesPLS<-function(X, Y, nComp=1, nrepeat=1000, sparsity='X', plot=TRUE, ask=FALSE, stepX=1, stepY=1){
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
  
  MaxVarX <- ncol(X)
  MaxVarY <- ncol(Y)
  
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



#################### Sparse PLS on X ####################

########## Models with 1 component

##### Calibration 

set.seed(1)
res=CalibratesPLS(X, Y, nComp = 1, nrepeat=1, sparsity='X')
## smallest MSEP around 0.98, with 16 variables in X

PlotCalibsPLS(res, sparsity='X', ncomp=1) # calibration plot


##### Run the model

MysPLSX <- spls(X,Y,keepX=1,ncomp=1,mode='regression')
MysPLSX$loadings$X # variables with loadings set to 0 are not selected


########## Models with 2 components

##### Calibration 

set.seed(1)
res=CalibratesPLS(X, Y, nComp = 2, nrepeat=1, sparsity='X')
# you have to enter manually the number of variables you want in the 1st component
# which is the number of variables leading to the smallest MSEP
# (It is done manually because if there are several low values, we put the smallest number of variables)
# In this instance, I put 16 for the 1st component and 16 again for the 2nd component


##### Run the model

MysPLSX2 <- spls(X,Y,keepX=c(1,6),ncomp=2,mode='regression')
MysPLSX2$loadings$X # one column by component


########## Sequential calibration (component by component)

##### 1st component 

set.seed(1)
res1=CalibratesPLSperComp(X, Y, comp = 1, nrepeat=1, sparsity='X') 

PlotCalibsPLS(res1, sparsity='X') 

SelectedX=4 # visual calibration


##### 2nd component 

res2=CalibratesPLSperComp(X, Y, comp = 2, SelectedX = SelectedX, nrepeat=1, sparsity='X', plot=FALSE, ask=FALSE)

res=rbind(res1,res2)
PlotCalibsPLS(res, sparsity='X', ncomp = 2) 

SelectedX=c(SelectedX, 5) # visual calibration


##### 3rd component

res3=CalibratesPLSperComp(X, Y, comp = 3, SelectedX = SelectedX, nrepeat=1, sparsity='X', plot=FALSE, ask=FALSE)

res=rbind(res,res3)
PlotCalibsPLS(res, sparsity='X', ncomp = 3) 

SelectedX=c(SelectedX, 1) # visual calibration



MysPLSX2 <- spls(X,Y,keepX=SelectedX,ncomp=3,mode='regression')
MysPLSX2$loadings$X # one column by component


#################### Sparse PLS on Y ####################

##### Calibration

set.seed(1)
res=CalibratesPLS(X, Y, nComp = 1, nrepeat=1, sparsity='Y')
## I selected 6


##### Run the model

MysPLSY <- spls(X,Y,keepY=6,ncomp=1,mode='regression')
MysPLSY$loadings$Y # we are interested in Y loadings


#################### Sparse PLS on X and Y ####################

##### Calibration

set.seed(1)
res=CalibratesPLS(X, Y, nComp = 1, nrepeat=1, sparsity='XY')
## This time the minimal MSEP values is assessed on a grid of 
## number of variables in X and Y
res[which.min(res$MSEP),]
## In this example, I chose 19 variables in 
## and 5 variables in Y

##### Run the model

MysPLSXY <- spls(X,Y,keepX=19, keepY=5, ncomp=1,mode='regression')
## We are interested in loadings for both X and Y:
MysPLSXY$loadings$X # in this particular exemple, all variables were selected (see calibration)
MysPLSXY$loadings$Y 






