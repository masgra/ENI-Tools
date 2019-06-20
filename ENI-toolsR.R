
####################################################################################
# -- General 
nc.path <- readLines( file("localpath.txt","r") ,n=1)
path.store <- paste(nc.path, "/Export/OTU-",sep="")

####################################################################################
# -- Preprocessing
source("preprocessing.R")


# load extracted files
otu <- read.csv(paste(path.store,"table-all.txt", sep=""), header = T, sep = "\t")

## 3) OTU rank based selection:
<<<<<<< HEAD
max.rank = 1000
=======
max.rank = 100
>>>>>>> e52fc565fdb368c565fe7787eeae46a9cbe9019f
min.rank = 1
otu.rank <- count_ranked_range(otu, min.lim = min.rank, max.lim = max.rank)
nOTU = max.rank

## 4) Subset separation and exclusion of rare species from OTU tabel 
Reactor = list(R.R = c("S3","S22","S27"), 
               R.C = c("S1","S4","S6","S8","S10","S12","S14","S16","S17","S19","S20","S23","S25"),
               R.D = c("S2","S5","S7","S9","S11","S13","S15","S18","S21","S24","S26")
)

# subset separation: reactor
otu.sub <- lapply(Reactor, function(x) otu[,which(colnames(otu) %in% x)])

# exclusion or rare species from subset 
# Note: **first exclusion, then separation by reactors would have been smarter
otu.sub <- exclude_rare(otu.sub, names.stay = otu.rank)


####################################################################################
# -- SPICE-EASI
library(SpiecEasi)
library(igraph)
library(Matrix)

# set number of replicates in StaRS
<<<<<<< HEAD
stars.replicates = 5
=======
stars.replicates = 20
>>>>>>> e52fc565fdb368c565fe7787eeae46a9cbe9019f

# exclude reference reactor
otu.sub$R.R <- NULL
SE <- list()

# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns
ptm <- proc.time()
# set parameters for glasso
<<<<<<< HEAD
param.glasso <- list("RC"=list(nlambda=50,"lambda"=1e-1),
                     "RD"=list(nlambda=50,"lambda"=1e-1))
=======
param.glasso <- list("RC"=list(nlambda=99,"lambda"=5e-1),
                     "RD"=list(nlambda=99,"lambda"=5e-1))
>>>>>>> e52fc565fdb368c565fe7787eeae46a9cbe9019f
# runn SPEC-EASI with glasso                   
SE$glasso <- mapply( function(x,param) {
  list(spiec.easi(t(x), nlambda= param$nlambda, method = "glasso", verbose = T, 
                  pulsar.params=list(rep.num=stars.replicates, ncores=1), 
                  lambda.min.ratio= param$lambda))}, 
  otu.sub, param.glasso)

<<<<<<< HEAD
# set parameters fro MB
param.mb <- list("RC"=list(nlambda=100,"lambda"=.5e-1),
                 "RD"=list(nlambda=100,"lambda"=.5e-1))
=======
a# set parameters fro MB
param.mb <- list("RC"=list(nlambda=99,"lambda"=5e-1),
                 "RD"=list(nlambda=99,"lambda"=5e-1))
>>>>>>> e52fc565fdb368c565fe7787eeae46a9cbe9019f
# runn SPEC-EASI with MB                   
SE$mb <- mapply( function(x,param) {
  list(spiec.easi(t(x), nlambda= param$nlambda, method = "mb", verbose = T, 
                  pulsar.params=list(rep.num=20, ncores=1), 
                  lambda.min.ratio= param$lambda))}, 
  otu.sub, param.mb)
proc.time() - ptm

# quanity check
# optimization: stability should be close to target stability 0.05. Increase nlambda 
# or lowering lambda.min.rat to get better target stability. However, if the network 
# gets empty (getOptInd <= 1), lambda.min.rate must be increased.  
lapply(SE$glasso, function(x) {list("0pt. lambda" = getOptLambda(x),
                                    "n of opt. lambda" = getOptInd(x),
                                    "stability" = getStability(x))})
lapply(SE$mb, function(x) {list("0pt. lambda" = getOptLambda(x),
                                "n of opt. lambda" = getOptInd(x),
                                "stability" = getStability(x))})

# sparCC: inital model estimation 
# Open:: 
# -- p-value calculation and network correction (not tested yet)
# -- add association sign (mutualExclusion, copresence) to igraph edges 


# Define (arbitrary) threshold for SparCC correlation matrix for the graph
sparcc.th= 0

SE$sparcc  <- lapply(otu.sub, function(x) {
  res <- sparcc(t(x), iter = 20, inner_iter = 30, th = 0.1)
  ## set values below threshold to zero
  refit <- abs(res$Cor) >= sparcc.th
  diag(refit) <- 0
  list(est = res, refit = refit) })

# SprarCC p-value estimation. 
# Multiply with refit matrix "p.refit" to apply significance filtering to the data.
compute.pval= FALSE
sparcc.p.th= .05
bootstraps = 100
if (compute.pval){
bootstap_p = function(x,y) {
    p.val <- pval.sparccboot(sparccboot(t(y), R=bootstraps, ncpus = 4 ), sided="both")
    pval <- matrix(0, round(.5+sqrt(2*length(p.val$pvals)+.25)), round(.5+sqrt(2*length(p.val$pvals)+.25)))
    pval[lower.tri(pval, diag=FALSE)] <- p.val$pvals
    p.refit <- as.matrix((t(pval)<=sparcc.p.th)*x$refit)
    p.refit[is.na(p.refit)] <-0
    return(list("est"=x$est, "refit"=x$refit, "p.values"=t(pval), "p.refit" = p.refit)) }
SE$sparcc$R.C <- bootstap_p(x=SE$sparcc$R.C, y=otu.sub$R.C)
SE$sparcc$R.D <- bootstap_p(x=SE$sparcc$R.D, y=otu.sub$R.C)
}

## create irgaph objects 
ig <- list()
ig$glasso <- lapply(SE$glasso, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(otu.sub$R.C))))
ig$mb <- lapply(SE$mb, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(otu.sub$R.C))))
ig$sparcc <- lapply(SE$sparcc, function(x) adj2igraph(x$refit,vertex.attr=list("name"=rownames(otu.sub$R.C))))
if(compute.pval){ # make p-value adgustment, if 
  ig$sparcc <- lapply(SE$sparcc, function(x) adj2igraph(x$p.refit,vertex.attr=list("name"=rownames(otu.sub$R.C))))
}


generat_edge_label <- function(graph, data, names){
  edges <- E(graph)
  param <- list()
  for(i in 1:length(edges)){
    nodes <- ends(graph,edges[i])
    param <- append(param,data[which(names==nodes[1]), which(names==nodes[2])])
  }
  return(param)
}


b <- lapply(a, function(x){
  if(x>0){
    return("copresence")
  }else if (x<0){
      return("mutualExclusion")
  }else{
    return("none")
  }})



a<- generat_edge_label(graph=ig$glasso$R.C,data=SE.res$glasso$R.C$weight ,names=rownames(otu.sub$R.C)  )



# Visualize using igraph plotting

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(t(otu.sub$R.D), 1))+4.5
# am.coord <- layout.fruchterman.reingold(ig.glasso)
am.coord <- layout.fruchterman.reingold(ig$mb$R.C)
# am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(2,3))
plot(ig$glasso$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig$mb$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
plot(ig$sparcc$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")
plot(ig$glasso$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig$mb$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
plot(ig$sparcc$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")


# https://kateto.net/networks-r-igraph





## set weights for final graph (possible threshold setting)
SE.res <-list()
SE.res$glasso <- lapply(SE$glasso, function(x) { list("weight"= as.matrix(cov2cor(getOptCov(x)) * getRefit(x))) })
SE.res$mb <- lapply(SE$mb, function(x) { list("weight"= as.matrix(as.matrix(getOptBeta(x)) * getRefit(x))) })
SE.res$sparcc <- lapply(SE$sparcc, function(x){ list("weight"= as.matrix(as.matrix(x$est$Cor) * x$refit)) })
if(compute.pval){ # make p-value adgustment, if 
  SE.res$sparcc <- lapply(SE$sparcc, function(x){ list("weight"= as.matrix(as.matrix(x$est$Cor) * x$p.refit)) })
}

hist(SE.res$glasso$R.C$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)
hist(SE.res$glasso$R.D$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)
hist(SE.res$mb$R.C$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)
hist(SE.res$mb$R.D$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)
hist(SE.res$sparcc$R.C$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)
hist(SE.res$sparcc$R.D$weight[SE.res$glasso$R.C$weight != 0], main='', xlab='edge weights', nclass =100)


# -- weigth filtering 
source("tools.R")

n.ranks = 50
how = "rank1" # ranks absolute values

SE.res$glasso <- lapply(SE.res$glasso, function(x) { list(weight= x$weight, th.weight = weight_filtering(x$weight, p=n.ranks, how=how))})
SE.res$mb <- lapply(SE.res$mb, function(x) { list(weight= x$weight, th.weight = weight_filtering(x$weight, p=n.ranks, how=how))})
SE.res$sparcc <- lapply(SE.res$sparcc, function(x) { list(weight= x$weight, th.weight = weight_filtering(x$weight, p=n.ranks, how=how))})


##  trasfere the networks to cytoscape

 if(!"RCy3" %in% installed.packages()){
   install.packages("BiocManager")
   BiocManager::install("RCy3")
 }
library(RCy3)

# better export matrix and use aMatReader 

# open cytoscape befor running because the network will be opened in cytoscape 3.x
createNetworkFromIgraph(ig$glasso$R.C,"RC-glasso")


a <- lapply(SE.res$glasso, function(x) { n <- data.frame(x$th.weight, row.names = names)
                                         colnames(n) <- names    
                                         return(n)})


w <-  lapply(SE.res, function(y) {
  lapply(y ,function(x) { 
    n <- data.frame(x$th.weight, row.names = names)
    colnames(n) <-names
    return(n)})
})


mapply(function(x,alg) {
  mapply(function(y,reac) {
    write.csv(y, paste(nc.path,"/Export2/weigths-",reac,"-",alg,"-all.csv", sep=""), quote = FALSE)},x ,names(x) )} ,w,names(w))
