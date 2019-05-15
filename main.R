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


# OTU rank based selection: get OTU sorted by mean ranking 
range <- 1:1250
otu.rank <- apply(otu, 2, function(x) rank(x)-min(rank(x)))
otu.rank <- rownames(otu.rank[which( rownames(otu.rank) %in% names(sort(rowSums(otu.rank),decreasing=T))[range]),])


## subset selection 

Reactor = list(R.R = c("S3", "S22", "S27"), 
               R.C = c("S1","S4","S6","S8","S10","S12","S14","S16","S17","S19","S20","S23","S25"),
               R.D = c("S2","S5","S7","S9","S11","S13","S15","S18","S21","S24","S26")
)

# subset selection: reactor
otu.sub <- otu[,which(colnames(otu) %in% Reactor$R.D)]
# store column sum 
otu.sub.sum <- colSums(otu.sub)
# subset selection: components 
otu.sub <- otu.sub[otu.rank, ]
# restore original size 
if ("Not_assigned" %in% otu.rank){
  # join removed counts to not assigend 
  otu.sub["Not_assigned",] <- otu.sub["Not_assigned",] + otu.sub.sum - colSums(otu.sub)
}else(
  # add new feature to OTU.sub with removed counts
  otu.sub[nrow(otu.sub)+1,] = list("Not_assigned" = otu.sub.sum - colSums(otu.sub))
)



# correlation between 
library("vegan")


plot(clr(otu.sub[,9]), clr(otu.sub[,11])) # high vs low expressed sample R.D
abline(a = 0, b = 1, col = 'red')

plot(clr(otu.sub[,11]), clr(otu.sub[,10])) # two low expressed sample R.D
abline(a = 0, b = 1, col = 'red')
colSums(otu.sub)

plot(clr(otu.sub[,2]), clr(otu.sub[,9])) # high vs low expressed sample R.C
abline(a = 0, b = 1, col = 'red')

plot(clr(otu.sub[,2]), clr(otu.sub[,10])) # two low expressed sample R.C
abline(a = 0, b = 1, col = 'red')






# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns
ptm <- proc.time()
SE.glasso  <- spiec.easi(t(otu.sub), nlambda=100, method = "glasso", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-4)
proc.time() - ptm

ptm <- proc.time()
SE.mb      <- spiec.easi(t(otu.sub), nlambda=100, method = "mb", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-4)
proc.time() - ptm

SE.sparcc  <- sparcc(t(otu.sub))

getOptInd(SE.glasso)
getStability(SE.glasso) # default target stability is 0.05 
getStability(SE.mb) # default target stability is 0.05 
sum(getRefit(SE.glasso))/2



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


