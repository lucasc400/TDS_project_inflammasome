rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)
library(parallel)

metab = readRDS(paste0(path,'metab_denoised_definite.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised.rds'))
#qt_score = readRDS(paste0(path,'proteins_qt_score.rds'))

# Functions
create_beta_matrix <- function(){
  M = matrix(0,dim(prot)[2],dim(metab)[2])
  colnames(M) = colnames(metab)
  rownames(M) = colnames(prot)
  return(M)
}

sel_pro_lasso.1se = create_beta_matrix()
dim(sel_pro_lasso.1se)

# Stability analysis
LassoSub = function(k = 1, Xdata, Ydata, family = "gaussian",
                    penalty.factor = NULL) {
  if (is.null(penalty.factor)) {
    penalty.factor = rep(1, ncol(Xdata))
  }
  set.seed(k)
  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata))
  Xsub = Xdata[s, ]
  Ysub = Ydata[s]
  model.sub = cv.glmnet(x = Xsub, y = Ysub, alpha = 1,
                        family = family, penalty.factor = penalty.factor)
  coef.sub = coef(model.sub, s = "lambda.1se")[-1]
  return(coef.sub)
}

LassoSub_loop = function(prot_vector, Xdata=metab){
  lasso.stab = sapply(1:100, FUN=LassoSub, Xdata=metab, Ydata=prot_vector)
  # IDEALLY - THIS OUTPUTS PROBABILITY. HOWEVER - INDEXING IS PROBLEM
  return()
}

niter = 100 # not used

# Parallelisation
nchunks = 20
ids = as.character(cut(1:ncol(prot),breaks=nchunks,labels=1:nchunks))

no_cores = detectCore
cl <- makeCluster(no_cores)
clusterExport(cl, c('LassoSub','LassoSub_loop','sel_pro_lasso.1se','metab','prot'))
clusterEvalQ(cl, library(glmnet))

t0 = Sys.time()

# PARALLEL CODE - INCOMPLETE
lasso.stab = parSapply(cl=cl, 1:nchunks, FUN=function(k){
  prot_chunk = prot[,ids==k]
  return(apply(prot_chunk, 2, FUN=LassoSub_loop))
})


sapply(1:20, FUN=function(k){
  prot_chunk = prot[,ids==k]
  as.numeric(apply(prot_chunk,2,mean))
})
V

# ORIGINAL CODE
for (i in 1:2){
  Y = prot[ ,i]
  lasso.stab = sapply(1:niter, FUN = LassoSub, Xdata = metab, Ydata = Y)
  lasso.prop = apply(lasso.stab, 1, FUN = function(x) {
    sum(x != 0)/length(x)
  })
  sel_pro_lasso.1se[i, ] = lasso.prop
}


stopCluster(cl)

t1 = Sys.time()
print(t1-t0)

#lasso.prop = apply(lasso.stab, 1, FUN = function(x) {sum(x != 0)/length(x)})


saveRDS(lasso.stab, paste0(path,'Lucas/Lasso_stability.rds'))
