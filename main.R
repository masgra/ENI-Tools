# isntall required packages
# install.packages("igraph") # network analysis tool
# install.packages("ggraph") # network analysis tool
# install.packages("vegan")  # pairwise dissimilarity 
# install.packages("MCL")    # cluster detection

# Install SpiecEasi package
# install.packages("devtools")
# library(devtools)
# install_github("zdk123/SpiecEasi", force = TRUE)

library(SpiecEasi)
library(igraph)
library(Matrix)

source("preprocessing.R")

nc.path <- readLines( file("localpath.txt","r") ,n=1)

## extract from biome file
library(phyloseq)

path.store <- paste(nc.path, "/Export/OTU-",sep="")
path.load <- paste(nc.path, "/Export/Comparison-EGGNOG.biom",sep="")

refresh.files = F
# refresh tabler files
if (refresh.files){
  otu <- biome2table(path.load, sum.dublicate.samples= TRUE, rename=TRUE, sort=TRUE)
  otu.t <- table2otu.table(otu$otu, otu$taxa)
  # write OTU table to file
  write.table(otu.t ,file = paste(path.store,"table-all.txt", sep=""),
              quote = T, eol='\n', na='NA', row.names = T, sep='\t')
  
  otu.q <- table2QUIIME.table(otu$otu, otu$taxa)
  # write OTU QUIIME table to file (for CoNet)
  write.table(otu.q,file = paste(path.store,"QUIIME-all.txt", sep=""),
              quote = F, sep='\t', eol='\n', na='NA', row.names = F)
  rm(otu, otu.t, otu.q)
}


## load extracted files
otu <- read.csv(paste(path.store,"table-all.txt", sep=""), header = T, sep = "\t")
otu$description <- NULL

## get signal-to-noise ratio
otu.sn <- sn_calc(otu, nrep=100, f.bootsp=0.01)





## subset selection 



























# set data path
path.data = readLines( file("localpath.txt","r") ,n=1)


## 1) Loading data into Phyloseq

# load feature counts
ds.all = read.csv(paste(path.data, "MT_Data-122018/readcounts/summary_L3.csv", sep=""), 
                  sep=";", header=T, row.names=1, check.names=F)
# remove unclassified reads
ds.all<-ds.all[-which(rownames(ds.all) == "Not"),]

# get sample 

Reactor = list(R.R = c("S-3", "S-22", "S-27"), 
               R.C = c("S-1","S-4","S-6","S-8","S-10","S-12","S-14","S-16","S-17","S-19","S-20","S-23","S-25"),
               R.D = c("S-2","S-5","S-7","S-9","S-11","S-13","S-15","S-18","S-21","S-24","S-26")
               )

# quick and dirty subset selection 
ds <- ds.all[,which(colnames(ds.all) %in% Reactor$R.C)]
otu.sub <- ds[which( rownames(ds) %in% names(sort(rowSums(ds),decreasing=T))[1:500]),]


# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns
ptm <- proc.time()
SE.glasso  <- spiec.easi(t(otu.sub), nlambda=10, method = "glasso", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-)
proc.time() - ptm
ptm <- proc.time()
SE.mb      <- spiec.easi(t(otu.sub), nlambda=100, method = "mb", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-4)
proc.time() - ptm
SE.sparcc  <- sparcc(t(otu.sub))

getOptInd(SE.mb)
getStability(SE.mb) # default target stability is 0.05 

sum(getRefit(SE.mb))/2



## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- ((abs(SE.sparcc$Cor) >= 0.3) & (abs(SE.sparcc$Cor) <= 0.99))
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

ig.glasso     <- adj2igraph(getRefit(SE.glasso))
ig.mb         <- adj2igraph(getRefit(SE.mb))
ig.sparcc     <- adj2igraph(sparcc.graph)


# Visualize using igraph plotting

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(otu.sub, 1))+4.5
am.coord <- layout.fruchterman.reingold(ig.glasso)
am.coord <- layout.fruchterman.reingold(ig.sparcc)
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.glasso, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")






secor  <- cov2cor(getOptCov(SE.glasso))
sebeta <- symBeta(getOptBeta(SE.mb), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(SE.glasso), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')


