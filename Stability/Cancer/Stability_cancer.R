rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
cov = readRDS(paste0(path,'covars_587.rds'))

metab_case = metab[!is.na(cov$subtype),]
metab_control = metab[is.na(cov$subtype),]

prot_case = prot[!is.na(cov$subtype),]
prot_control = prot[is.na(cov$subtype),]

# Set up output matrix
sel_pro_case = matrix(0, dim(prot)[2], dim(metab_case)[2])
colnames(sel_pro_case) = colnames(metab_case)
rownames(sel_pro_case) = colnames(prot)
sel_pro_control = matrix(0, dim(prot)[2], dim(metab_control)[2])
colnames(sel_pro_control) = colnames(metab_control)
rownames(sel_pro_control) = colnames(prot)

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
  # Case
  Y_case = prot_case[ ,i]
  lasso.stab.case = sapply(1:niter, FUN = LassoSub, Xdata = metab_case, Ydata = Y_case)
  lasso.prop.case = apply(lasso.stab.case, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_case[i, ] = lasso.prop.case
  # Control
  Y_control = prot_control[, i]
  lasso.stab.control = sapply(1:niter, FUN = LassoSub, Xdata = metab_control, Ydata = Y_control)
  lasso.prop.control = apply(lasso.stab.control, 1, FUN = function(x) {sum(x != 0)/length(x)})
  sel_pro_control[i, ] = lasso.prop.control
}

t1 = Sys.time()
print(t1-t0)



saveRDS(sel_pro_case, '../results/Lucas/Stability/Cancer/sel_pro_case.rds')
saveRDS(sel_pro_control, '../results/Lucas/Stability/Cancer/sel_pro_control.rds')
