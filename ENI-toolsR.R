
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
max.rank = 20
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
stars.replicates = 50

# exclude reference reactor
otu.sub$R.R <- NULL
SE <- list()

# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns
ptm <- proc.time()
# set parameters for glasso
param.glasso <- list("RC"=list(nlambda=100,"lambda"=1e-1),
                     "RD"=list(nlambda=50,"lambda"=1e-1))
# runn SPEC-EASI with glasso                   
SE$glasso <- mapply( function(x,param) {
  list(spiec.easi(t(x), nlambda= param$nlambda, method = "glasso", verbose = T, 
                  pulsar.params=list(rep.num=stars.replicates, ncores=1), 
                  lambda.min.ratio= param$lambda))}, 
  otu.sub, param.glasso)

# set parameters fro MB
param.mb <- list("RC"=list(nlambda=50,"lambda"=1e-1),
                 "RD"=list(nlambda=50,"lambda"=1e-1))
# runn SPEC-EASI with MB                   
SE$mb <- mapply( function(x,param) {
  list(spiec.easi(t(x), nlambda= param$nlambda, method = "mb", verbose = T, 
                  pulsar.params=list(rep.num=stars.replicates, ncores=1), 
                  lambda.min.ratio= param$lambda))}, 
  otu.sub, param.mb)
proc.time() - ptm

# quanity check
# optimization: stability should be close to target stability 0.05. Increase nlambda 
# or lambda.min.rat to get better target stability. However, if the network 
# gets empty (getOptInd <= 1), lambda.min.rate must be increased.  
lapply(SE$glasso, function(x) {list("0pt. lambda" = getOptLambda(x),
                                    "n of opt. lambda" = getOptInd(x),
                                    "stability" = getStability(x))})
lapply(SE$mb, function(x) {list("0pt. lambda" = getOptLambda(x),
                                "n of opt. lambda" = getOptInd(x),
                                "stability" = getStability(x))})



Sparcc.res  <- lapply(otu.sub, function(x) sparcc(t(x), iter = 20, inner_iter = 30, th = 0.25))
# p-value calcualtion for sparcc
Sparcc.pval <- lapply(otu.sub, function(x) pval.sparccboot(sparccboot(t(x), R =1000, ncpus=2), sided = "both"))
# save sparcc p-value file 
save(Sparcc.pval, file=paste(nc.path ,"Export2/sparccPVAL-N1000-01.RData",sep=""))

vertix <- as.list(name=rownames(otu.sub$R.C)); 
names(vertix) <- rownames(otu.sub$R.C)

## Define (arbitrary) threshold for SparCC correlation matrix for the graph
SE$sparcc <- lapply(Sparcc.res, function(x) (abs(x$Cor) >= 0.9) & (abs(x$Cor) < 1.0))
lapply(SE$sparcc, function(x) diag(x) <-0)
SE$sparcc <- lapply(SE$sparcc, function(x) Matrix(x, sparse=TRUE))

# Save the SE object
save(SE, file=paste(nc.path ,"Export2/igraph-N1000-01.RData",sep=""))

## create irgaph objects 
ig <- list()
ig$glasso <- lapply(SE$glasso, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(otu.sub$R.C))))
ig$mb <- lapply(SE$mb, function(x) adj2igraph(getRefit(x),vertex.attr=list("name"=rownames(otu.sub$R.C))))
ig$sparcc <- lapply(SE$sparcc, function(x) adj2igraph(x,vertex.attr=list("name"=rownames(otu.sub$R.C))))



# Visualize using igraph plotting

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(t(otu.sub$R.D), 1))+4.5
# am.coord <- layout.fruchterman.reingold(ig.glasso)
am.coord <- layout.fruchterman.reingold(ig$sparcc$R.C)
# am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(2,3))
plot(ig$glasso$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig$mb$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
plot(ig$sparcc$R.C, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")
plot(ig$glasso$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig$mb$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
plot(ig$sparcc$R.D, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")


# https://kateto.net/networks-r-igraph








secor  <- cov2cor(getOptCov(SE.glasso))
sebeta <- symBeta(getOptBeta(SE.mb), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(SE.glasso), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')


##  trasfere the networks to cytoscape

# if(!"RCy3" %in% installed.packages()){
#   install.packages("BiocManager")
#   BiocManager::install("RCy3")
# }
library(RCy3)

# open cytoscape befor running because the network will be opened in cytoscape 3.x
createNetworkFromIgraph(ig$glasso$R.C,"RC-glasso")
createNetworkFromIgraph(ig$glasso$R.D,"RD-glasso")
createNetworkFromIgraph(ig$mb$R.C,"RC-MB")
createNetworkFromIgraph(ig$mb$R.D,"RD-MB")
createNetworkFromIgraph(ig$sparcc$R.C,"RC-sparcc")
createNetworkFromIgraph(ig$sparcc$R.D,"RD-sparcc")