# isntall required packages
# install.packages("igraph") # network analysis tool
# install.packages("ggraph") # network analysis tool
# install.packages("vegan")  # pairwise dissimilarity 
# install.packages("MCL")    # cluster detection

# Install SpiecEasi package
# install.packages("devtools")
# library(devtools)
# install_github("zdk123/SpiecEasi", force = TRUE)

library(igraph)
library(Matrix)

####################################################################################
# -- General 
nc.path <- readLines( file("localpath.txt","r") ,n=1)




####################################################################################
# -- Preprocessing

source("preprocessing.R")


## 1) extract from biome file
library(phyloseq)

# shoul biom files be extracted again  
refresh.files = F

path.store <- paste(nc.path, "/Export/OTU-",sep="")
path.load <- paste(nc.path, "/Export/Comparison-EGGNOG.biom",sep="")

# refresh and load OTU tabler files

if (refresh.files){
  otu <- biome2table(path.load, sum.dublicate.samples= TRUE, rename=TRUE, sort=TRUE)
  otu.t <- table2otu.table(otu$otu, otu$taxa)
  # write OTU table to file
  write.table(otu.t ,file = paste(path.store,"table-all.txt", sep=""),
              quote = T, eol='\n', na='NA', row.names = T, sep='\t')
  
  otu.q <- table2QUIIME.table(otu$otu, otu$taxa)
  
  #replace 'not assigend' - line
  otu.q[which(otu.q$`#OTU ID` == "-2"),ncol(otu.q)]
  
  # write OTU QUIIME table to file (for CoNet)
  write.table(otu.q,file = paste(path.store,"QUIIME-all.txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  otu <- otu.t
  rm(otu.t, otu.q)
}else{
  # load extracted files
  otu <- read.csv(paste(path.store,"table-all.txt", sep=""), header = T, sep = "\t")
}


## 2) Get signal-to-noise ratio (we do not use it because it is computationaly heavy)
# otu.sn <- sn_calc(otu, nrep=100, f.bootsp=0.01)


## 3) OTU rank based selection:
max.rank = 1000
min.rank = 1
otu.rank <- count_ranked_range(otu, min.lim = min.rank, max.lim = max.rank)
nOTU = max.rank

## 4) Subset separation and exclusion of rare species from OTU tabel 
Reactor = list(R.R = c("S3","S22","S27"), 
               R.C = c("S1","S4","S6","S8","S10","S12","S14","S16","S17","S19","S20","S23","S25"),
               R.D = c("S2","S5","S7","S9","S11","S13","S15","S18","S21","S24","S26")
)

# subset separation: reactor
otu.sub <- list()
otu.sub$R.C <- otu[,which(colnames(otu) %in% Reactor$R.C)]
otu.sub$R.D <- otu[,which(colnames(otu) %in% Reactor$R.D)]
otu.sub$R.R <- otu[,which(colnames(otu) %in% Reactor$R.R)]

# exclusion or rare species from subset 
# Note: **first exclusion, then separation by reactors would have been smarter
otu.sub <- exclude_rare(otu.sub, names.stay = otu.rank)


## 5) Subset separation and exclusion of rare species from QUIIME tabel

reset.file = F
if (reset.file){
  otu.q <- read.csv(paste(path.store,"QUIIME-all.txt", sep=""), header = T, sep = "\t", 
                    check.names = TRUE, row.names = 1)
  otu.q <- exclude_rare_QUIIME(otu.q, names.stay = otu.rank, resetID = T)
  
  descripts <- c('taxonomy', '#OTU ID')
  
  otu.q.sub <- list()
  otu.q.sub$R.C <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.C, descripts))]
  otu.q.sub$R.D <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.D, descripts))]
  otu.q.sub$R.R <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.R, descripts))]
  
  otu.q.sub<-lapply(otu.q.sub,function(x) cbind( (x[,1:(ncol(x)-1)]/colSums(x[,1:(ncol(x)-1)])),x[,ncol(x)]))
  
  # write OTU QUIIME table to file (for CoNet)
  write.table(cbind('#OTU ID' = rownames(otu.q.sub$R.C), otu.q.sub$R.C) ,file = paste(path.store,"QUIIME-RC-rar-abs",nOTU,".txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  write.table(cbind('#OTU ID' = rownames(otu.q.sub$R.D), otu.q.sub$R.D) ,file = paste(path.store,"QUIIME-RD-rar",nOTU,".txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  write.table(cbind('#OTU ID' = rownames(otu.q.sub$R.R), otu.q.sub$R.R) ,file = paste(path.store,"QUIIME-RR-rar",nOTU,".txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  
  write.table(cbind('#OTU ID' = rownames(otu.q), otu.q) ,file = paste(path.store,"QUIIME-rar",nOTU,".txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  
  rm(otu.q, otu.q.sub, descripts)
}



####################################################################################
# -- general inspection of DS

l <- nrow(otu.sub$R.C)-1

# correlation between 
library("SpiecEasi")
library("vegan")

# start vs end 
par(mfrow=c(1,3))
plot(clr(otu.sub$R.C$S1[1:l]), clr(otu.sub$R.C$S25[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')
plot(clr(otu.sub$R.R$S3[1:l]), clr(otu.sub$R.R$S27[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')
plot(clr(otu.sub$R.D$S2[1:l]), clr(otu.sub$R.D$S26[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')


# high vs. low library size 
par(mfrow=c(1,3))
plot(clr(otu.sub$R.C$S25[1:l]), clr(otu.sub$R.C$S10[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')
plot(clr(otu.sub$R.C$S25[1:l]), clr(otu.sub$R.C$S4[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')
plot(clr(otu.sub$R.D$S9[1:l]), clr(otu.sub$R.D$S11[1:l])) # high vs low expressed sample R.D
  abline(a = 0, b = 1, col = 'red')

  
  
# get inverse simpson indices 
vegan::diversity(as.matrix(otu[,2:ncol(otu)]), index = "invsimpson", MARGIN = 2)
vegan::diversity(as.matrix(otu.sub$R.C), index = "invsimpson", MARGIN = 2) 
vegan::diversity(as.matrix(otu.sub$R.D), index = "invsimpson", MARGIN = 2) 
vegan::diversity(as.matrix(otu.sub$R.R), index = "invsimpson", MARGIN = 2) 


rm(l)


####################################################################################
# -- plot frequencies of metrics
metric <- list()
pcount <- .5

# calculate relative abundancies
rab.sub <- lapply(otu.sub, function(x) {
  # add small pseudo-count
  x[x<1.0] <- pcount
  res <- sweep(x,2,colSums(x),"/")
  return(res)})
  
# Pearson Correlation
metric$pearson <- lapply(rab.sub, function(x) {
  res <- cor(t(x),method = "pearson")
  diag(res)=NA
  # remove "not asignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not_assigned")),])
  res$"Not_assigned" <- NULL
  return(as.matrix(res)) })

# Spearman correlation
metric$spearman <- lapply(rab.sub, function(x) {
  res <- cor(t(x),method = "spearman")
  diag(res)=NA
  # remove "not asignede" (trash component) component
  res <-  data.frame(res[which(!(rownames(res) == "Not_assigned")),])
  res$"Not_assigned" <- NULL
  return(as.matrix(res))})

## scaled variance log-ratio
# metric$VLR <- lapply(rab.sub, function(x){
#   res <- 1-exp(-sqrt(balance::vlr(t(x), alpha = min(x)/2)))
#   diag(res)=NA
#   # remove "not asignede" (trash component) component
#   res <-  as.data.frame(res[which(!(rownames(res) == "Not_assigned")),])
#   res$"Not_assigned" <- NULL
#   return(as.matrix(res)) })

# Kullback-Leibler dissimilarity
metric$KL <- lapply(rab.sub, function(x){
  res <- flexmix::KLdiv(t(x),eps=1e-8)
  res <- res + t(res)
  diag(res)=NA
  # remove "not asignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not_assigned")),])
  res$"Not_assigned" <- NULL
  return(as.matrix(res)) })

#Bray Curtis dissimilarity
metric$BC<-lapply(rab.sub,function(x) {
  res <- sweep(x,1,rowSums(x),"/") # set rowsum to 1 
  res <- as.matrix(vegan::vegdist(res, method="bray", binary=F, diag=T, upper=T, na.rm = T))
  diag(res)=NA
  # remove "not asignede" (trash component) component
  res <-  as.data.frame(res[which(!(rownames(res) == "Not_assigned")),])
  res$"Not_assigned" <- NULL
  return(as.matrix(res)) })

# get metric distributions 
metric.dist <- rapply(metric, function(x) c(x[lower.tri(x,diag = F)]), how = 'list')
rm(metric)

# get upper and lower quantiles for nvertex = 100 
n.vert <- (nrow(rab.sub$R.C)-1)*(nrow(rab.sub$R.C)-2)/2
n.th <- 499 # = n=1000!
metric.th <- c((n.vert-n.th)/n.vert, 1-(n.vert-n.th)/n.vert)

# compute quantils: 
th.quntils <- rapply(metric.dist, function(x) as.matrix(quantile( x , probs = metric.th)) , how = 'list')

# write threshold file
source("tools.R")
methods.names= c('correl_pearson', 'correl_spearman', #'sim_varlogratio', 
                 'dist_kullbackleibler', 'dist_bray')
write_th_file(th.quntils,methods.names,paste(nc.path,'/CoNet-IN/th_file-union-', sep=''))


# plot frequency histograms
library("ggplot2")
#Pearson
ggplot()+ 
  stat_density( aes( x=metric.dist$pearson$R.C, colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=metric.dist$pearson$R.D, colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Pearson Coefficient Dencity', x="pairwise Pearson Coefficient", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-1,1))

#Spearman
ggplot()+ 
  stat_density( aes( x=metric.dist$spearman$R.C, colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=metric.dist$spearman$R.D, colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Spearman Coefficient Dencity', x="pairwise Spearman Coefficient", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-1,1))

# #VLR
# ggplot()+ 
#   stat_density( aes( x=c(na.omit(c(metric$VLR$R.C))), colour="R.C"), na.rm = T, size=1, alpha=.2)+
#   stat_density( aes( x=c(na.omit(c(metric$VLR$R.D))), colour="R.D"), na.rm = T, size=1, alpha=.2)+
#   scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
#   labs(title ='Pairwise Variance of log-ratio Dencity', x="pairwise variance of log-ratios", y="densiy" )+
#   theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
#   scale_x_continuous(limits=c(0.,1.0))

#KL
ggplot()+ 
  stat_density( aes( x=c(na.omit(c(metric$KL$R.C))), colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=c(na.omit(c(metric$KL$R.D))), colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pairwise Kullback-Leibler dissimilarity Density', x="pairwise Kullback-Leibler dissimilarity", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,8))

#BC
ggplot()+ 
  stat_density( aes( x=c(na.omit(c(metric$BC$R.C[lower.tri(metric$BC$R.C,diag = F)]))), colour="R.C"), na.rm = T, size=1, alpha=.2)+
  stat_density( aes( x=c(na.omit(c(metric$BC$R.D[lower.tri(metric$BC$R.D,diag = F)]))), colour="R.D"), na.rm = T, size=1, alpha=.2)+
  scale_color_manual("Reactro", values = c("red", "blue"), labels = c(bquote(~R[C]), bquote(~R[D])))+
  labs(title ='Pariwise Bray Curtis dissimilarity Density', x="pairwise Bray Curtis dissimilarity", y="densiy" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(0,1))



## compuet ranks: low values get smallest ranks
metric.r <- list()
metric.r$pearson <- lapply(metric.dist$pearson, function(x) length(x) - rank( x ,ties.method = "average" , na.last = T)) 
metric.r$spearman <- lapply(metric.dist$spearman, function(x) length(x) - rank(x,ties.method = "average" , na.last = T)) 
# metric.r$VLR <- lapply(metric.dist$VLR, function(x) rank(x , ties.method = "average" , na.last = T)) 
metric.r$KL <- lapply(metric.dist$KL, function(x) rank(x,ties.method = "average" , na.last = T))
metric.r$BC <- lapply(metric.dist$BC, function(x) rank(x,ties.method = "average" , na.last = T))

# calculate union rank
metric.r$union <- NULL
n.seq <- seq_along(metric.r)
for (i in names(metric.r$pearson)) {
  metric.r$union[[i]]<- as.matrix(rowSums(as.data.frame( lapply(n.seq, 
                                          function(j) unlist(metric.r[[j]][i])))))
}
rm(n.seq)


# get upper and lower quantiles for nvertex = 100 
n.vert <- (nrow(rab.sub$R.C)-1)*(nrow(rab.sub$R.C)-2)/2
n.th <- 499 # = n=1000!
metric.th <- c((n.vert-n.th)/n.vert, 1-(n.vert-n.th)/n.vert)


# get postionons that exceed thresholds
metric.exceed.pos <- lapply(metric.r$union, function(x){
  q <- quantile(x, probs = metric.th)
  res <- list(high.rank = which(x >= q[1]), 
              low.rank  = which(x <= q[2]))
  return(res)
})

# calculate thresholds from min/max values of exceeding positions

bb <- list()

bb$pearson <- lapply(names(metric.dist$pearson), function(x){
  c( min( metric.dist$pearson[[x]][metric.exceed.pos[[x]]$low.rank], na.rm = T), 
     max( metric.dist$pearson[[x]][metric.exceed.pos[[x]]$high.rank], na.rm = T))}) 
names(bb$pearson) <- names(metric.dist$pearson)

bb$spearman <- lapply(names(metric.dist$spearman), function(x){
  c( min( metric.dist$spearman[[x]][metric.exceed.pos[[x]]$low.rank], na.rm = T), 
     max( metric.dist$spearman[[x]][metric.exceed.pos[[x]]$high.rank], na.rm = T))}) 
names(bb$spearman) <- names(metric.dist$spearman)

# bb$VLR <- lapply(names(metric.dist$VLR), function(x){
#   c( min( metric.dist$VLR[[x]][metric.exceed.pos[[x]]$high.rank], na.rm = T), 
#      max( metric.dist$VLR[[x]][metric.exceed.pos[[x]]$low.rank], na.rm = T))}) 
# names(bb$VLR) <- names(metric.dist$VLR)

bb$KL <- lapply(names(metric.dist$KL), function(x){
  c( min( metric.dist$KL[[x]][metric.exceed.pos[[x]]$high.rank], na.rm = T), 
     max( metric.dist$KL[[x]][metric.exceed.pos[[x]]$low.rank], na.rm = T))}) 
names(bb$KL) <- names(metric.dist$KL)

bb$BC <- lapply(names(metric.dist$BC), function(x){
  c( min( metric.dist$BC[[x]][metric.exceed.pos[[x]]$high.rank], na.rm = T), 
     max( metric.dist$BC[[x]][metric.exceed.pos[[x]]$low.rank], na.rm = T))}) 
names(bb$BC) <- names(metric.dist$BC)


# write threshold file
source("tools.R")
methods.names= c('correl_pearson', 'correl_spearman', #'sim_varlogratio', 
                 'dist_kullbackleibler', 'dist_bray')
write_th_file(bb,methods.names,paste(nc.path,'/CoNet-IN/th_file-intersection-', sep=''))



rm(bb,metric.dist,metric.exceed.pos,metric.r,th.quntils, i, methods.names,metric.th, min.rank, n.vert)

  ####################################################################################
  # -- SPICE-EASI
library("SpiecEasi")


ds <- otu.sub$R.C
  
  
# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns

# optimization: if 

SE.glasso <- list()
SE.mb <- list()
SE.sparcc <- list()

ptm <- proc.time()
SE.glasso$RC  <- spiec.easi(t(otu.sub$R.C), nlambda=50, method = "glasso", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-1)
SE.glasso$RD  <- spiec.easi(t(otu.sub$R.D), nlambda=50, method = "glasso", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-1)
proc.time() - ptm
# quanity check 
getOptInd(SE.glasso$RC)# get index of optimum lambda: close to 1 is good for comutational speed, exactly 1 is bad. 
getOptInd(SE.glasso$RD)
getStability(SE.glasso$Rc) # achieved stability: target stability is 0.05 
getStability(SE.glasso$RD)
# optimization: 

ptm <- proc.time()
SE.mb$RC      <- spiec.easi(t(otu.sub$R.C), nlambda=50, method = "mb", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-1)
SE.mb$RD      <- spiec.easi(t(otu.sub$R.D), nlambda=50, method = "mb", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-1)
proc.time() - ptm
proc.time() - ptm
# quality check
getOptInd(SE.mb$RC) 
getOptInd(SE.mb$RD) 
getStability(SE.mb$RC)
getStability(SE.mb$RD)


SE.sparcc$RC  <- sparcc(otu.sub$R.C, iter = 20, inner_iter = 20, th = 0.1)
SE.sparcc$RD  <- sparcc(otu.sub$R.D, iter = 20, inner_iter = 20, th = 0.1)

vertix <- as.list(rownames(ds)); names(vertix) <- rownames(ds)

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph.RC <- ((abs(SE.sparcc$RC$Cor) >= 0.3) & (abs(SE.sparcc$RC$Cor) < 1.0))
sparcc.graph.RD <- ((abs(SE.sparcc$RD$Cor) >= 0.3) & (abs(SE.sparcc$RD$Cor) < 1.0))
diag(sparcc.graph.RC) <- 0
diag(sparcc.graph.RD) <- 0

sparcc.graph.RC <- Matrix(sparcc.graph.RC, sparse=TRUE)
sparcc.graph.RD <- Matrix(sparcc.graph.RD, sparse=TRUE)

ig.glasso <- list(); ig.mb <- list(); ig.sparcc <- list() 

ig.glasso$RC     <- adj2igraph(getRefit(SE.glasso$RC), vertex.attr= vertix)
ig.mb$RC         <- adj2igraph(getRefit(SE.mb$RC))
ig.sparcc$RC     <- adj2igraph(sparcc.graph.RC)

ig.glasso$RD     <- adj2igraph(getRefit(SE.glasso$RD), vertex.attr= vertix)
ig.mb$RD         <- adj2igraph(getRefit(SE.mb$RD))
ig.sparcc$RD     <- adj2igraph(sparcc.graph,RD)


# Visualize using igraph plotting

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(otu.sub$R.C, 1))+4.5
# am.coord <- layout.fruchterman.reingold(ig.glasso)
am.coord <- layout.fruchterman.reingold(ig.sparcc$R.C)
# am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(2,3))
  plot(ig.glasso$RC, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  plot(ig.mb$RC, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
  plot(ig.sparcc$RC, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")
  plot(ig.glasso$RD, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  plot(ig.mb$RD, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
  plot(ig.sparcc$RD, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")


# https://kateto.net/networks-r-igraph


plot(1,2)






secor  <- cov2cor(getOptCov(SE.glasso))
sebeta <- symBeta(getOptBeta(SE.mb), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(SE.glasso), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')


