# This script generates the threshold files for CoNet. 
# currently, therer are four different metrics from CoNet implemented in this script. 
# All metrics measure pairwise (univariate) associations. The thresholds limits are set
# based on the number of pairwise assuciations to consider. 
# Two different types of threshold files are implemented:
#
# union: each metric result is ranked. Limits of the upper and the lower quantiles 
# are chosen as thresholds. Threshold files can be generated in parallel for multiple 
# datasets (experiments). (This procedure controls the number of assocaitions per metric)
# 
# intersection: each metric result is ranked. Ranks of all metrics are sumed up. 
# Qunatiles of the sumed ranks are computed and by that, the set of components. The
# thresholds are set to the highes/lowes value of the corresponding metrics values of
# this set of component. (This procedure controls the absolute number of associations for 
# which all metrics agree)
# 



####################################################################################
# -- General

source("tools.R")
# library("SpiecEasi")
library("vegan")

## set paths (create output directories)
# input directory
path.in <- './Data'
# output directory
dir.create(file.path(".", "Results"),showWarnings = FALSE) #does not overwrite an existing directory
path.out <- './Results'
dir.create(file.path(path.out, "CoNet"),showWarnings = FALSE) #does not overwrite an existing directory

## load dataset
load(file=paste(path.in,"/MT-OTU-Reactors-1000.RData", sep=""))
ds <- otu.r.sub
rm(otu.r.sub)
# remove reference reactor
ds$R.R <- NULL



####################################################################################
# -- inspection of DS

# compute inverse simpson indices 
vegan::diversity(as.matrix(ds$R.C), index = "invsimpson", MARGIN = 2) 
vegan::diversity(as.matrix(ds$R.D), index = "invsimpson", MARGIN = 2) 


####################################################################################
# -- plot frequency distribution of metrics
metric <- list()
pcount <- 0 # pseudo count (not necessary)

# calculate relative abundancies, add pseudo counts (only to zeros)
ds.r <- lapply(ds, function(x) {
  # add small pseudo-count
  x[x<1.0] <- pcount
  res <- sweep(x,2,colSums(x),"/")
  return(res)})

# compute Pearson Correlation
metric$correl_pearson <- lapply(ds.r, function(x) {
  res <- cor(t(x),method = "pearson")
  diag(res)=NA
  # remove "not assignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not")),])
  res$"Not" <- NULL
  return(as.matrix(res)) })

# compute Spearman correlation
metric$correl_spearman <- lapply(ds.r, function(x) {
  res <- cor(t(x),method = "spearman")
  diag(res)=NA
  # remove "not assignede" (trash component) component
  res <-  data.frame(res[which(!(rownames(res) == "Not")),])
  res$"Not" <- NULL
  return(as.matrix(res))})


## scaled variance log-ratio
## this metric has been omitted due to a bug in CoNet
# metric$VLR <- lapply(ds.r, function(x){
#   res <- 1-exp(-sqrt(balance::vlr(t(x), alpha = min(x)/2)))
#   diag(res)=NA
#   # remove "not assignede" (trash component) component
#   res <-  as.data.frame(res[which(!(rownames(res) == "Not")),])
#   res$"Not" <- NULL
#   return(as.matrix(res)) })

# compute Kullback-Leibler dissimilarity
metric$dist_kullbackleibler <- lapply(ds.r, function(x){
  res <- flexmix::KLdiv(t(x),eps=1e-8)
  res <- res + t(res)
  diag(res)=NA
  # remove "not assignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not")),])
  res$"Not" <- NULL
  return(as.matrix(res)) })

# compute Bray Curtis dissimilarity
metric$dist_bray<-lapply(ds.r,function(x) {
  res <- sweep(x,1,rowSums(x),"/") # set rowsum to 1 
  res <- as.matrix(vegan::vegdist(res, method="bray", binary=F, diag=T, upper=T, na.rm = T))
  diag(res)=NA
  # remove "not assignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not")),])
  res$"Not" <- NULL
  return(as.matrix(res)) })

# compute distribution vectors for each metric (and experiment)
metric.dist <- rapply(metric, function(x) 
  c(x[lower.tri(x,diag = F)]), how = 'list')
rm(metric)

## plot frequency histograms
library("ggplot2")
#Pearson
ggplot()+ 
  stat_density( aes( x=metric.dist$correl_pearson$R.C, colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=metric.dist$correl_pearson$R.D, colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Pearson Coefficient Dencity', x="pairwise Pearson Coefficient", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-1,1))

#Spearman
ggplot()+ 
  stat_density( aes( x=metric.dist$correl_spearman$R.C, colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=metric.dist$correl_spearman$R.D, colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Spearman Coefficient Dencity', x="pairwise Spearman Coefficient", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-1,1))

#KL
ggplot()+ 
  stat_density( aes( x=c(na.omit(c(metric.dist$dist_kullbackleibler$R.C))), colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=c(na.omit(c(metric.dist$dist_kullbackleibler$R.D))), colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Kullback-Leibler dissimilarity Density', x="pairwise Kullback-Leibler dissimilarity", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,8))

#BC
ggplot()+ 
  stat_density( aes( x=c(na.omit(c(metric.dist$dist_bray$R.C[lower.tri(metric.dist$dist_bray$R.C,diag = F)]))), colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=c(na.omit(c(metric.dist$dist_bray$R.D[lower.tri(metric.dist$dist_bray$R.D,diag = F)]))), colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pariwise Bray-Curtis dissimilarity Density', x="pairwise Bray-Curtis dissimilarity", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,1))



####################################################################################
# --  Compute threshold and write threshold files

## rank the metric distributions: low values -> samll ranks
metric.r <- list()
metric.r <- lapply(metric.dist, function(x){
  lapply(x, function(y){
    rank( y ,ties.method = "average" , na.last = T)}) 
  })

# swap the ranks of pearson and spearman: high valures -> samll ranks
metric.r$correl_pearson <-lapply(metric.r$correl_pearson, function(x) length(x)-x)
metric.r$correl_spearman <-lapply(metric.r$correl_spearman, function(x) length(x)-x)

# calculate sum rank vector(s)
rank.sum <- list()
for (i in names(metric.r[[1]])){
  print(i)
  rank.sum[[i]] <- rowSums(as.data.frame(lapply(metric.r, function(x) x[i])))
}

# calculate the upper and lower edge limit, relative to the total number of possible 
# edges. 
# i.e.: (500 upper, 500 lower for n.th=1000 with total edges= 498501 --> [99.9%,  00.1%])
# total number of edges per metric
n.th <- 1000 
n.vert <- (nrow(ds.r$R.C)-1)*(nrow(ds.r$R.C)-2)/2
# relative edge quantil
metric.th <- c((n.vert-(n.th/2))/n.vert, 1-(n.vert-(n.th/2))/n.vert) 



## write union threshold files

# compute metric quantils 
th.quantils.union <- rapply(metric.dist, function(x) quantile( x , probs = metric.th) , how = 'list')
# write threshold file for CoNet
writeThresholdFile(th.quantils.union,path.out=paste(path.out,'/CoNet/th_file-union-',n.th,"-", sep=''))



## write intersection threshold files

# compute list of associations from upper/lower rank sum quantiles 
rank.sum.pos <- lapply(rank.sum, function(x){
  q <- quantile(x, probs = metric.th)
  res <- list(mut.ex = which(x >= q[1]), 
              cooc  = which(x <= q[2]))
  return(res)
})

# calculate thresholds: min/max of metric values from those associations within
# the upper/lower rank sum quantiles. Note that pearson & spearman again must be swaped!
th.quantils.intersect <- list()

th.quantils.intersect$correl_pearson <- mapply(function(x,y){
   a <- data.frame(c(min( x[y$cooc], na.rm = T), 
         max( x[y$mut.ex], na.rm = T)))
   colnames(a) <- names(x)
   return(a)
  }, x= metric.dist$correl_pearson, y= rank.sum.pos)

th.quantils.intersect$correl_spearman <- mapply(function(x,y){
  a <- data.frame(c(min( x[y$cooc], na.rm = T), 
               max( x[y$mut.ex], na.rm = T)))
  colnames(a) <- names(x)
  return(a)
}, x= metric.dist$correl_spearman, y= rank.sum.pos)

th.quantils.intersect$dist_kullbackleibler <- mapply(function(x,y){
  a <- data.frame(c(min( x[y$mut.ex], na.rm = T),
             max( x[y$cooc], na.rm = T)))
  colnames(a) <- names(x)
  return(a)
}, x= metric.dist$dist_kullbackleibler, y= rank.sum.pos)

th.quantils.intersect$dist_bray <- mapply(function(x,y){
  a <- data.frame(c(min( x[y$mut.ex], na.rm = T),
          max( x[y$cooc], na.rm = T)))
  colnames(a) <- names(x)
  return(a)
}, x= metric.dist$dist_bray, y= rank.sum.pos)


# write threshold file
writeThresholdFile(th.quantils.intersect,path.out=paste(path.out,'/CoNet/th_file-intersect-',n.th,"-", sep=''))



