# PreProcessing. This script:
#
# - reduces the number of components in the count table to the top xxx while keeping the 
#   absolute number of counts (keep the closure).
# - splits the count table into subsets (experiments)
# - replaces zeros with a pseudo counts


####################################################################################
# -- General

source("tools.R")

## set paths
# input directory
path.in <- './Data'

# read count table
otu <- read.csv2(file=paste(path.in,"MT-OTU-all.txt",sep = "/"), 
                  check.names = FALSE, sep="\t", header = TRUE, row.names = 1)
# set pseudo count 
p.count = 0



## OTU rank based selection settinges
max.rank = 1000
min.rank = 1
# get name vectro of selected components
otu.rank <- SelectRankedComponents(otu, min.lim = min.rank, max.lim = max.rank)
# get number of components
nOTU = max.rank-min.rank+1



## Excluding rare components without omiting the counts (preserves the closure)
otu.r.all <- excludeComponents(otu, names.stay=otu.rank, trash.node.name="Not_assigned")
# replace zeros with a pseudo count
otu.r.all[otu.r.all<=0] <- p.count



## Dividing the count table into experiment subsets 
Reactor = list(R.R = c("3","22","27"), 
               R.C = c("1","4","6","8","10","12","14","16","17","19","20","25"),
               R.D = c("2","5","7","9","11","13","15","18","21","24","26"))
# split dataset into reactor subsets
otu.r.sub <- lapply(Reactor, function(x){
  otu.r.all[,which(colnames(otu.r.all) %in% x)]
})
# write count tables for each reactro to a separate file
write.table(otu.r.sub$R.C ,file = paste(path.in,"/MT-OTU-RC-",nOTU,".txt", sep=""),
            quote = F, sep='\t', eol='\n', na='NA', row.names = T, col.names = T)
write.table(otu.r.sub$R.D ,file = paste(path.in,"/MT-OTU-RD-",nOTU,".txt", sep=""),
            quote = F, sep='\t', eol='\n', na='NA', row.names = T, col.names = T)
# store the subset count table object
save(otu.r.sub, file=paste(path.in,"/MT-OTU-Reactors-",nOTU,".RData", sep=""))





####################################################################################
## generate CoNet input table: QIIME table of condensated count table for each experiment (reactor)

# get ortholog structure
ortholog <- read.csv2(file=paste(path.in,"MT-ortholog-all.txt",sep = "/"), 
                  check.names = FALSE, sep="\t", header = TRUE, row.names = 1)
# rearange ortholog structure to QIIME fromat 
ortholog <- data.frame( taxonomy= paste(ortholog$L0,ortholog$L1,ortholog$L2,ortholog$L3, sep = "; "), row.names= rownames(ortholog))

# merge filtered count tables and ortholog column
ds <- list()
ds$R.C <- merge(x = cbind("#OTU"=rownames(otu.r.sub$R.C),otu.r.sub$R.C), y = ortholog , by.x = "#OTU" , by.y = "row.names", all.x = TRUE)
ds$R.D <- merge(x = cbind("#OTU"=rownames(otu.r.sub$R.D),otu.r.sub$R.D), y = ortholog , by.x = "#OTU" , by.y = "row.names", all.x = TRUE)

# write QIIME tables for each reactro to a separate file
write.table(ds$R.C ,file = paste(path.in,"/MT-QIIME-RC-",nOTU,".txt", sep=""),
            quote = F, sep='\t', eol='\n', na='NA', row.names = F, col.names = T)
write.table(ds$R.D ,file = paste(path.in,"/MT-QIIME-RD-",nOTU,".txt", sep=""),
            quote = F, sep='\t', eol='\n', na='NA', row.names = F, col.names = T)



