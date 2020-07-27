rm(list=ls())
path = '/rdsgpfs/general/project/hda_tds/live/GROUP_2/used/'
setwd(path)

metab = readRDS(paste0(path,'metab_denoised_definite.rds'))
metab = as.matrix(metab)
metab = t(metab)
prot = readRDS(paste0(path,'proteins_denoised.rds'))
qt_score = readRDS(paste0(path,'proteins_qt_score.rds'))

suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(factoextra))
library(dplyr)

# Check the scales of data

colmean=colMeans(metab)
colmean=as.data.frame(colmean)
ggplot(aes(x=colmean),data=colmean)+geom_density() # indeed mean is 0

colsd=apply(metab,1,sd)
colsd=as.data.frame(colsd)
summary(colsd) #min=1113, median=11183, max=6205988

metab_df = as.data.frame(metab)
metab_df %>%
  summarise(mean, sd)

# PCA and clustering tendancy
#fviz_pca_ind(prcomp(metab),geom='point')

#res = get_clust_tendency(metab, n=nrow(metab)-1, graph=FALSE)
#res$hopkins_stat # 0.8645117

# K-means clustering
# elbow plot
fviz_nbclust(metab,kmeans,method='wss') # suggests three clusters as best
# silhouette plot
fviz_nbclust(metab, kmeans, method='silhouette') # also suggests three clusters

metab_scaled = scale(metab)
# elbow plot
fviz_nbclust(metab_scaled,kmeans,method='wss') # no clear cut-point
# silhouette plot
fviz_nbclust(metab_scaled, kmeans, method='silhouette') # suggests two clusters

# Is it the case that the clusters found by k-means are based on sizes of the values?

metab_kmeans = kmeans(metab,centers = 3, nstart=20)
metab_kmeans$cluster
