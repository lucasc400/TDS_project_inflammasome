rm(list=ls())
#path=dirname(rstudioapi::getActiveDocumentContext()$path)
#path=str_remove(path,'/Scripts/.*')
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)
library(glmnet)

metab = readRDS(paste0(path,'metab_1776.rds'))
metab = as.matrix(metab)
prot = readRDS(paste0(path,'proteins_denoised_new.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score_new.rds')

# Univariate LASSO proteins regressed with metabolites
create_beta_matrix <- function(){
  M = matrix(0,dim(prot)[2],dim(metab)[2])
  colnames(M) = colnames(metab)
  rownames(M) = colnames(prot)
  return(M)
}

#beta_lasso_lambda.min = create_beta_matrix()
beta_lasso_lambda.1se = create_beta_matrix()
#beta_en_lambda.min = create_beta_matrix()
#beta_en_lambda.1se = create_beta_matrix()


t0 = Sys.time()
for (i in 1:dim(prot)[2]){
  y = prot[ ,i]
  model_itr = cv.glmnet(metab, y, alpha=1, family='gaussian')
#  beta_lasso_lambda.min[i, ] = coef(model_itr, s=model_itr$lambda.min)[-1]
  beta_lasso_lambda.1se[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
  print(paste('lambda.min =', model_itr$lambda.min, 'lambda.1se =', model_itr$lambda.1se))
}
t1 =  Sys.time()
print(t1-t0)

#apply(beta_lasso_lambda.1se,1,function(x) sum(x!=0))
metab_sel_count = apply(beta_lasso_lambda.1se, 2, function(x){sum(x!=0)})
metab_sel_count = sort(metab_sel_count,decreasing=TRUE)

write_rds(beta_lasso_lambda.1se, paste0(path, '../results/Lucas/metab_beta_lasso.rds'))
write_rds(metab_sel_count, paste0(path, '../results/Lucas/metab_sel_count.rds'))


# LASSO stratified by smoking ------------------------------------------------
cov = readRDS(paste0(path,'covars_587.rds'))
head(cov)
cov$case =ifelse(is.na(cov$subtype),0,1)
cov$smoking = ifelse(cov$smoking_status=='Never',0,1)

metab_never = metab[cov$smoking==0,]
metab_ever = metab[cov$smoking==1,]

prot_never = prot[cov$smoking==0,]
prot_ever = prot[cov$smoking==1,]

beta_never = matrix(0, 92, dim(metab_never)[2])
colnames(beta_never) = colnames(metab_never)
rownames(beta_never) = colnames(prot_never)

beta_ever = matrix(0, 92, dim(metab_ever)[2])
colnames(beta_ever) = colnames(metab_ever)
rownames(beta_ever) = colnames(prot_ever)


t0 = Sys.time()
for (i in 1:92){
  y_never = prot_never[,i]
  y_ever = prot_ever[,i]

  model_itr = cv.glmnet(metab_never, y_never, alpha=1, family='gaussian')
  beta_never[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
  
  model_itr = cv.glmnet(metab_ever, y_ever, alpha=1, family='gaussian')
  beta_ever[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
}
t1 = Sys.time()
print(t1-t0)

beta_smoking <- list('ever'=beta_ever,'never'=beta_never)
saveRDS(beta_smoking, paste0(path,'../results/Lucas/beta_smoking.rds'))

# LASSO stratified by case ------------------------------------------------
metab_control = metab[cov$case==0,]
metab_case = metab[cov$case==1,]

prot_control = prot[cov$case==0,]
prot_case = prot[cov$case==1,]

beta_control = matrix(0, 92, dim(metab_control)[2])
colnames(beta_control) = colnames(metab_control)
rownames(beta_control) = colnames(prot_control)

beta_case = matrix(0, 92, dim(metab_case)[2])
colnames(beta_case) = colnames(metab_case)
rownames(beta_case) = colnames(prot_case)


t0 = Sys.time()
for (i in 1:92){
  y_control = prot_control[,i]
  y_case = prot_case[,i]
  
  model_itr = cv.glmnet(metab_control, y_control, alpha=1, family='gaussian')
  beta_control[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
  
  model_itr = cv.glmnet(metab_case, y_case, alpha=1, family='gaussian')
  beta_case[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
}
t1 = Sys.time()
print(t1-t0)

# Output results
beta_cancer <- list('case'=beta_case,'control'=beta_control)
saveRDS(beta_cancer, paste0(path,'../results/Lucas/beta_cancer.rds'))

# Elastic net (decided not to be used) ------------------------------------

cvm.1se = function(alpha) {
  model = cv.glmnet(x=metab, y=Y, alpha=alpha)
  with(model, cvm[which.min(lambda-lambda.1se)])
}

cvm.min = function(alpha) {
  model = cv.glmnet(x=metab, y=Y, alpha=alpha)
  with(model, cvm[which.min(lambda-lambda.min)])
}



t0 = Sys.time()
for (i in 1:dim(prot)[2]){
  Y = prot[ ,i]
  alpha.opt = optimise(cvm.1se, c(0,1))
  model_itr = cv.glmnet(metab, y, alpha=alpha.opt$minimum, family='gaussian')
  beta_en_lambda.min[i, ] = coef(model_itr, s=model_itr$lambda.min)[-1]
  beta_en_lambda.1se[i, ] = coef(model_itr, s=model_itr$lambda.1se)[-1]
  print(paste('lambda.min =', model_itr$lambda.min, 'lambda.1se =', model_itr$lambda.1se))
}
t1 =  Sys.time()

