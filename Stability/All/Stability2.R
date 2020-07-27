rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)
library(parallel)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score_new.rds'))

# Functions
create_beta_matrix <- function(){
  M = matrix(0,dim(prot)[2]/4,dim(metab)[2])
  colnames(M) = colnames(metab)
  return(M)
}

sel_pro_lasso.1se = create_beta_matrix()
rownames(sel_pro_lasso.1se) = colnames(prot)[24:46]
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

niter = 100

# Parallelisation
#cl <- makeCluster(no_cores)
#clusterExport(cl, c('LassoSub','sel_pro_lasso.1se','metab','prot'))

# prot - Y matrix, 92 columns

t0 = Sys.time()
for (i in 1:23){
  Y = prot[ ,i+23]
  lasso.stab = sapply(1:niter, FUN = LassoSub, Xdata = metab, Ydata = Y)
  lasso.prop = apply(lasso.stab, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_lasso.1se[i, ] = lasso.prop
}

t1 = Sys.time()
print(t1-t0)



saveRDS(sel_pro_lasso.1se, paste0(path,'../results/Lucas/Stability/All/sel_pro_2_new.rds'))
