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
max.rank = 250
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

reset.file = FALSE
if (reset.file){
  otu.q <- read.csv(paste(path.store,"QUIIME-all.txt", sep=""), header = T, sep = "\t", 
                    check.names = TRUE, row.names = 1)
  otu.q <- exclude_rare_QUIIME(otu.q, names.stay = otu.rank, resetID = T)
  
  descripts <- c('taxonomy', '#OTU ID')
  
  otu.q.sub <- list()
  otu.q.sub$R.C <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.C, descripts))]
  otu.q.sub$R.D <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.D, descripts))]
  otu.q.sub$R.R <- otu.q[ ,which(colnames(otu.q) %in% c(Reactor$R.R, descripts))]
  
  
  # write OTU QUIIME table to file (for CoNet)
  write.table(cbind('#OTU ID' = rownames(otu.q.sub$R.C), otu.q.sub$R.C) ,file = paste(path.store,"QUIIME-RC-rar",nOTU,".txt", sep=""),
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

  rm(l)

  ####################################################################################
  # -- SPICE-EASI
library("SpiecEasi")

ds <- otu.sub$R.C
  
  
# Spiec easy: non-normalized count OTU/data table with samples on rows and features/OTUs in columns

# optimization: if 


ptm <- proc.time()
SE.glasso  <- spiec.easi(t(ds), nlambda=50, method = "glasso", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-2)
proc.time() - ptm
# quanity check 
getOptInd(SE.glasso) # get index of optimum lambda: close to 1 is good for comutational speed, exactly 1 is bad. 
getStability(SE.glasso) # achieved stability: target stability is 0.05 
# optimization: 

ptm <- proc.time()
SE.mb      <- spiec.easi(t(ds), nlambda=50, method = "mb", verbose = T, pulsar.params=list(rep.num=10, ncores=1), lambda.min.ratio=1e-1)
proc.time() - ptm
# quality check
getOptInd(SE.mb) 
getStability(SE.mb)


SE.sparcc  <- sparcc(ds, iter = 20, inner_iter = 20, th = 0.1)

vertix <- as.list(rownames(ds)); names(vertix) <- rownames(ds)

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- ((abs(SE.sparcc$Cor) >= 0.3) & (abs(SE.sparcc$Cor) < 1.0))
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

ig.glasso     <- adj2igraph(getRefit(SE.glasso), vertex.attr= vertix)
ig.mb         <- adj2igraph(getRefit(SE.mb))
ig.sparcc     <- adj2igraph(sparcc.graph)


# Visualize using igraph plotting

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(ds, 1))+4.5
am.coord <- layout.fruchterman.reingold(ig.glasso)
am.coord <- layout.fruchterman.reingold(ig.sparcc)
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
  plot(ig.glasso, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
  plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="mb")
  plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Sparcc")


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


