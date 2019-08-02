


####################################################################################
# -- General

source("tools.R")
library(SpiecEasi)
library(igraph)
library(Matrix)



## set paths (create output directories)
# input directory
path.in <- './Data'
# output directory
dir.create(file.path(".", "Results"),showWarnings = FALSE) #does not overwrite an existing directory
path.out <- './Results'
dir.create(file.path(path.out, "SpiecEasi"),showWarnings = FALSE) #does not overwrite an existing directory

## load dataset
load(file=paste(path.in,"/MT-OTU-Reactors-1000.RData", sep=""))
ds <- otu.r.sub
rm(otu.r.sub)
# remove reference reactor
ds$R.R <- NULL


####################################################################################
# -- SPICE-EASI: compute glasso, MB and sparCC

# set memory limit (in my case - Windows - oversizing is okay)
# memory.limit(16108*2)

# set number of StaRS replicates. This is the number of bootstraps that StaRS picks to 
# evaluate lambda. Increase this number if possible. 
stars.replicates = 20

# result list
SE <- list()

# flag to indicate whether to refreshing results, if results already exist
refresh = FALSE

## run and store glasso results - if not already exist or if refresh is forced
if(!file.exists(paste(path.out, "/SpiecEasi/glasso.RData", sep=""))|refresh ){
  # set parameters for glasso
  param.glasso <- list("RC"=list(nlambda=20,"lambda"=5e-1),
                       "RD"=list(nlambda=20,"lambda"=7.5e-1))
  # runn SPEC-EASI with glasso 
  SE$glasso <- mapply( function(x,param) {gc()
    list(spiec.easi(t(x), nlambda= param$nlambda, method = "glasso", verbose = T, 
                    pulsar.params=list(rep.num=stars.replicates, ncores=1), 
                    lambda.min.ratio= param$lambda))}, 
    ds, param.glasso)
  # store glasso results
  glasso <- SE$glasso
  save( glasso , file = paste(path.out, "/SpiecEasi/glasso.RData", sep=""))
  rm(glasso)
}else{
  load(paste(path.out, "/SpiecEasi/glasso.RData", sep=""))
  SE$glasso <- glasso
  rm(glasso)
}

# run and store MB results - if not already exist or if refresh is forced
if(!file.exists(paste(path.out, "/SpiecEasi/MB.RData", sep=""))|refresh ){

  # set parameters fro MB
  param.mb <- list("RC"=list(nlambda=20,"lambda"=5e-2),
                   "RD"=list(nlambda=20,"lambda"=5e-2))
  # runn SPEC-EASI with MB                   
  SE$mb <- mapply( function(x,param) {gc()
    list(spiec.easi(t(x), nlambda= param$nlambda, method = "mb", verbose = T, 
                    pulsar.params=list(rep.num=stars.replicates, ncores=1), 
                    lambda.min.ratio= param$lambda))}, 
    ds, param.mb)
  # store mb results
  mb <- SE$mb
  save( mb , file = paste(path.out, "/SpiecEasi/MB.RData", sep=""))
  rm(mb)
}else{
  load(paste(path.out, "/SpiecEasi/MB.RData", sep=""))
  SE$mb <- mb
  rm(mb)
}



## quanity check
# optimization: stability should be close to target stability 0.05. Increase nlambda 
# or lambda.min.rat to get better target stability. However, if the network 
# gets empty (getOptInd <= 1), lambda.min.rate must be increased.  
lapply(SE$glasso, function(x) {list("0pt. lambda" = getOptLambda(x),
                                    "n of opt. lambda" = getOptInd(x),
                                    "stability" = getStability(x))})
lapply(SE$mb, function(x) {list("0pt. lambda" = getOptLambda(x),
                                "n of opt. lambda" = getOptInd(x),
                                "stability" = getStability(x))})

# sparCC: inital model estimation 
# --> p-value calculation and network correction is !!untested!!


## run and store sparCC results - if not already exist or if refresh is forced
if(!file.exists(paste(path.out, "/SpiecEasi/sparCC.RData", sep=""))|refresh ){
  
  # arbitrary threshold set to ignree very uncorrelated edges 
  sparcc.th= 0 # ignore threshold
  SE$sparcc  <- lapply(ds, function(x) {
    res <- sparcc(t(x), iter = 20, inner_iter = 20, th = 0.1)
    # set values below threshold to zero
    refit <- abs(res$Cor) >= sparcc.th
    # set diagonal to zero
    diag(refit) <- 0
    list(est = res, refit = refit) })
  # store spraCC results
  sparCC <- SE$sparcc
  save( sparCC , file = paste(path.out, "/SpiecEasi/sparCC.RData", sep=""))
}else{
  load(paste(path.out, "/SpiecEasi/sparCC.RData", sep=""))
  SE$sparcc <- sparCC
}

# !!UNTESTED!! sparCC p-value estimation. 
# Multiply with refit matrix "p.refit" to apply significance filtering to the data.
compute.pval= FALSE
if ((!file.exists(paste(path.out, "/SpiecEasi/sparCC-pval.RData", sep=""))|refresh) & compute.pval){
  sparcc.p.th= .001
  bootstraps = 1000
  bootstap_p = function(x,y) {
    p.val <- pval.sparccboot(sparccboot(t(y), R=bootstraps, ncpus = 4 ), sided="both")
    pval <- matrix(0, round(.5+sqrt(2*length(p.val$pvals)+.25)), round(.5+sqrt(2*length(p.val$pvals)+.25)))
    pval[lower.tri(pval, diag=FALSE)] <- p.val$pvals
    p.refit <- as.matrix((t(pval)<=sparcc.p.th)*x$refit)
    p.refit[is.na(p.refit)] <-0
    return(list("est"=x$est, "refit"=x$refit, "p.values"=t(pval), "p.refit" = p.refit)) }
  sparcc$R.C <- bootstap_p(x=SE$sparcc$R.C, y=ds$R.C)
  sparcc$R.D <- bootstap_p(x=SE$sparcc$R.D, y=ds$R.C)
  # OPEN::
  # - store results of p-value calculation separately
  # - create full p-value matrix (not just upper triangular)
}



####################################################################################
# -- visualize initial networks 

# # create irgaph objects
# ig <- list()
# ig$glasso <- lapply(SE$glasso, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(ds$R.C))))
# ig$mb <- lapply(SE$mb, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(ds$R.C))))
# ig$sparcc <- lapply(SE$sparcc, function(x) adj2igraph(x$refit,vertex.attr=list("name"=rownames(ds$R.C))))
# if(compute.pval){ # make p-value adgustment, if
#   ig$sparcc <- lapply(SE$sparcc, function(x) adj2igraph(x$p.refit,vertex.attr=list("name"=rownames(ds$R.C))))
# }
# 
# # Visualize using igraph
# 
# ## set size of vertex proportional to clr-mean
# vsize$R.C    <- rowMeans(clr(t(ds$R.C), 1))+4.5
# vsize$R.D    <- rowMeans(clr(t(ds$R.D), 1))+4.5
# am.coord <- layout.fruchterman.reingold(ig$sparcc$R.C)
# 
# par(mfrow=c(2,3))
# plot(ig$glasso$R.C, layout=am.coord, vertex.size=vsize$R.C, vertex.label=NA, main="glasso-RC")
# plot(ig$mb$R.C, layout=am.coord, vertex.size=vsize$R.C, vertex.label=NA, main="mb-RC")
# plot(ig$sparcc$R.C, layout=am.coord, vertex.size=vsize$R.C, vertex.label=NA, main="Sparcc-RC")
# plot(ig$glasso$R.D, layout=am.coord, vertex.size=vsize$R.D, vertex.label=NA, main="glasso-RD")
# plot(ig$mb$R.D, layout=am.coord, vertex.size=vsize$R.D, vertex.label=NA, main="mb-RD")
# plot(ig$sparcc$R.D, layout=am.coord, vertex.size=vsize$R.D, vertex.label=NA, main="Sparcc-RD")





####################################################################################
# -- inspect and define network weights


## calculate weights of associations
SE.res <-list()
# glasso: uses correlation as weights 
SE.res$glasso <- lapply(SE$glasso, function(x) { list("weight"= as.matrix(cov2cor(getOptCov(x)) * getRefit(x))) })
# MB: uses beta as weights 
SE.res$mb <- lapply(SE$mb, function(x) { list("weight"= as.matrix(as.matrix(symBeta(getOptBeta(x),mode='ave')) * getRefit(x))) })
# sparCC: uses correlation as weights 
SE.res$sparcc <- lapply(SE$sparcc, function(x){ list("weight"= as.matrix(as.matrix(x$est$Cor) * x$refit)) })
# for sparCC apply p-value adjusted refit, if exists 
# (sets all associations with p-values larger than "sparcc.p.th" to zero)
if(!is.null(SE$sparcc[1]$p.refit)){
  SE.res$sparcc <- lapply(SE$sparcc, function(x){ list("weight"= as.matrix(as.matrix(x$est$Cor) * x$p.refit)) })
}


## create histograms of weights distribution. Zeros are omitted.
par(mfrow=c(2,3))
hist(SE.res$sparcc$R.C$weight[abs(SE.res$sparcc$R.C$weight) > 0.0], main='sparCC-RC', xlab='edge weights', nclass =100)
hist(SE.res$glasso$R.C$weight[abs(SE.res$glasso$R.C$weight) > 0.0], main='gLasso-RC', xlab='edge weights', nclass =100)
hist(SE.res$mb$R.C$weight[abs(SE.res$mb$R.C$weight) > 0.0], main='MB-RC', xlab='edge weights', nclass =100)

hist(SE.res$sparcc$R.D$weight[abs(SE.res$sparcc$R.D$weight) > 0.0], main='sparCC-RD', xlab='edge weights', nclass =100)
hist(SE.res$glasso$R.D$weight[abs(SE.res$glasso$R.D$weight) > 0.0], main="gLasso-RD", xlab='edge weights', nclass =100)
hist(SE.res$mb$R.D$weight[abs(SE.res$mb$R.D$weight) > 0.0], main='MB-RD', xlab='edge weights', nclass =100)




####################################################################################
# -- set and apply weight filtering 

# set limits for number of edges 
n.ranks = list("R.C" = 1000, 
               "R.D" = 1000)
# identifyer for weightFiltering()
#  "rank1": omit all except for highest absolut association weights. 
#  "rank2": omit all except for highest and lowest association weights (half/half). 
#  "percentil": same as "rank2", but handles input as percentils instead of number of associations.  
how = "rank1" 

# compute weight filtered adjacency matrix
SE.res <- lapply(SE.res, function(z) { 
  mapply(function(x,y) {
    return( list("weight"= x$weight, "th.weight" = weightFiltering(x$weight, p=y, how=how)))}, 
    z, n.ranks, SIMPLIFY = FALSE, USE.NAMES = TRUE)})




####################################################################################
# -- export networks to edge tables 

source("tools.R")

# this function converts the interaction type, based on weights (definition equivalent to CoNet)
getInteraction = function(x , inter=c("copresence","none","mutualExclusion")){
  r <- c()
  r[which(x>0)] <- inter[1]
  r[which(x==0)] <- inter[2]
  r[which(x<0)] <- inter[3]
  return(r)
}

# calculate relative abundance median per component and reactor (I added this as a vertex feature)
vertex.features <- lapply(ds, function(x) list("median.abd"=apply( sweep(x,2,colSums(x),`/`) ,1,FUN = median)))


## add phylo-structure to export
# this requires the structuere conversion of the tree export from MEGAN into a table structure
# The script for this conversion also part of this project, but was implemented in python! 
phylo <- read.csv(paste(path.in, "/MT-ortholog-all.txt", sep=""), header = T, sep = "\t", stringsAsFactors = F)

# create list of featuere data frame
vertex.features <- lapply(vertex.features, function(x) cbind(phylo, "median.abd"= x$median.abd[match(rownames(phylo), names(x$median.abd))]))



# export associations to files (one file per algorithm and per experiment)
mapply(function(x,alg) {
  mapply(function(y,reac,z) {
    exp <- getCytoscapeEdgeDS(y, names=rownames(ds$R.C), expand.Unreported.Verteces = T, feature = z)
    write.table(exp, paste(path.out,"/SpiecEasi/weigths-",reac,"-",alg,".txt", sep=""),sep = "\t", quote = F, na = "",row.names = F, dec = ".")},
    x , names(x), vertex.features )} ,
  SE.res , names(SE.res) )







