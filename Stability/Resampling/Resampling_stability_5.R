# Resampling for stability
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score_new.rds'))


set.seed(5)
subset = sample(1:587, 293)
metab_subset = metab[subset,]
prot_subset = prot[subset,]

create_beta_matrix <- function(){
  M = matrix(0,dim(prot)[2],dim(metab)[2])
  colnames(M) = colnames(metab)
  return(M)
}

sel_pro_lasso.1se = create_beta_matrix()
rownames(sel_pro_lasso.1se) = colnames(prot_subset)

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

t0 = Sys.time()
for (i in 1:92){
  Y = prot_subset[ ,i]
  lasso.stab = sapply(1:niter, FUN = LassoSub, Xdata = metab_subset, Ydata = Y)
  lasso.prop = apply(lasso.stab, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_lasso.1se[i, ] = lasso.prop
}

t1 = Sys.time()
print(t1-t0)

print(sum(sel_pro_lasso.1se!=0))
saveRDS(sel_pro_lasso.1se, paste0(path,'../results/Lucas/Stability/Resampling/sel_pro_r5.rds'))