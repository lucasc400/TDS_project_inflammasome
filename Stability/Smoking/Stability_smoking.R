rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
cov = readRDS(paste0(path,'covars_587.rds'))

metab_ever = metab[cov$smoking_status!='Never',]
metab_never = metab[cov$smoking_status=='Never',]

prot_ever = prot[cov$smoking_status!='Never',]
prot_never = prot[cov$smoking_status=='Never',]

# Set up output matrix
sel_pro_ever = matrix(0, dim(prot)[2], dim(metab_ever)[2])
colnames(sel_pro_ever) = colnames(metab_ever)
rownames(sel_pro_ever) = colnames(prot)
sel_pro_never = matrix(0, dim(prot)[2], dim(metab_never)[2])
colnames(sel_pro_never) = colnames(metab_never)
rownames(sel_pro_never) = colnames(prot)

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
  # Ever
  Y_ever = prot_ever[ ,i]
  lasso.stab.ever = sapply(1:niter, FUN = LassoSub, Xdata = metab_ever, Ydata = Y_ever)
  lasso.prop.ever = apply(lasso.stab.ever, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_ever[i, ] = lasso.prop.ever
  # Never
  Y_never = prot_never[, i]
  lasso.stab.never = sapply(1:niter, FUN = LassoSub, Xdata = metab_never, Ydata = Y_never)
  lasso.prop.never = apply(lasso.stab.never, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_never[i, ] = lasso.prop.never
}

t1 = Sys.time()
print(t1-t0)



saveRDS(sel_pro_ever, '../results/Lucas/Stability/Smoking/sel_pro_ever.rds')
saveRDS(sel_pro_never, '../results/Lucas/Stability/Smoking/sel_pro_never.rds')
