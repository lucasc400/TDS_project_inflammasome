rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score_new.rds'))

# Univariate LASSO SIS regressed with metabolites

model_lasso = cv.glmnet(metab, qt_score, alpha=1, family='gaussian')
plot(model_lasso)

model_lasso$lambda.1se
length(which(coef(model_lasso, s=model_lasso$lambda.1se)!=0)) #6 metabolites

betas_sis = coef(model_lasso, s='lambda.1se')[-1]
names(betas_sis) = rownames(coef(model_lasso, s='lambda.1se'))[-1]

# unstable


# Stability ---------------------------------------------------------------

# Functions 
sel_pro_qt = matrix(0, 1, dim(metab)[2])
colnames(sel_pro_qt) = colnames(metab)

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
lasso.stab = sapply(1:niter, FUN = LassoSub, Xdata = metab, Ydata = qt_score)
lasso.prop = apply(lasso.stab, 1, FUN = function(x) {sum(x != 0)/length(x)})
sel_pro_qt[1, ] = lasso.prop

dim(sel_pro_qt)
t1 = Sys.time()
print(t1-t0)

sel_pro_qt_df = as.data.frame(sel_pro_qt)
colnames(sel_pro_qt_df) = colnames(metab)
sum(sel_pro_qt_df!=0)
sort(sel_pro_qt_df,decreasing=TRUE)[1:69]

saveRDS(sel_pro_qt_df,'../results/Lucas/Stability/qt_stability.rds')
