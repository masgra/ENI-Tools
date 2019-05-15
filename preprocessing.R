


# calculate signal to noise ratio
sn_calc <- function(x, nrep=10000, f.bootsp=0.01, addCount =0.25, f.scale=1, mc=TRUE) {
  
  x <- (round(otu*f.scale)+addCount)
  
  if (mc){
    library('parallel')
    
    # generate cluster 
    cl <- makeCluster(detectCores()-1); clusterEvalQ(cl,library(MASS))
    # Set a different seed on each member of the cluster (just in case)
    clusterSetRNGStream(cl)
    # set input (just in case)
    # clusterExport(cl,c("x", "f.bootsp", "nrep"))
    #... then parallel replicate...
    system.time({
      #  set input (this is one column of x)
      a <- lapply(x, function(x){
        # put objects in place that might be needed for the code
        #clusterExport(cl,c("x", "f.bootsp", "nrep"))
        parSapply(cl, 1:nrep, # build replicates 
                    function(i,...) {
                      tabulate(  # sum counts 
                        sample( 1:length(x), size= round(sum(x)*f.bootsp), replace=T, prob = x/sum(x) ) # resampling with replacement
                        )})
      }) 
    })
    #stop the cluster
    stopCluster(cl)
  }else{
    # single core version
    system.time({
      a <- lapply(x , 
                  function(x){  data.frame( replicate(nrep, tabulate( 
                     sample( 1:length(x), size= round(sum(x)*f.bootsp), replace=T, prob = x/sum(x) ))))}    )
     })
  }
  
    
  # get mean of replicates per OTU  
  otu.mean <- as.data.frame(lapply(a, function(x) rowMeans(x)))
  # get variance of replicates per OTU
  otu.var <- as.data.frame(lapply(a, function(x) apply(x, 1, FUN = var)))
  # build signal-to-noise table
  sn.ratio <- otu.mean/otu.var
  rownames(sn.ratio) <- rownames(otu)
  
  return (list("nrep" = nrep, 
               "f.bootsp" = f.bootsp, 
               "addCount" = addCount, 
               "f.scale" = f.scale,
               "replicates" = a,
               "otu.mean" = otu.mean,
               "otu.var" = otu.var, 
               "otu.snratio" = sn.ratio))
}


biome2table <- function(filepath, sum.dublicate.samples= TRUE, rename=TRUE, sort=TRUE){ 
  
  if ((rename==F)&&(sort==T)){
    warning('"sort=TRUE" is ignored because "rename=FALES"')
  }
  
  # load data form biom-format file
  data <- phyloseq::import_biom(filepath)
  # get OTU table as data frame
  otu <- data.frame(data@otu_table@.Data, check.names=F, stringsAsFactors = FALSE)
  # get taxa, levels separated as columns
  taxa <-data.frame(data@tax_table@.Data, stringsAsFactors = FALSE)
  if(sum.dublicate.samples){
    # sum samples with same ID
    otu <- sapply( unique( stringr::str_extract(colnames(otu), "ID-[0-9]+-") ), 
                   function(x) rowSums(otu[,grep(x, colnames(otu)), drop=FALSE]))
  }
  # rearange order of columns to increasing by sample number  
  if(rename){
    colnames(otu) <- stringr::str_extract(colnames(otu), "[0-9]+") # rename
    if (sort){
      otu <- otu[,order(colnames(otu))] # sort 
    }
    colnames(otu) <- paste('S',colnames(otu), sep = '') # rename
  }
  
  return (list("otu" = otu,
               "taxa" = taxa))
}


# create table from MEGAN6 biom file export of lowest taxa level 
table2otu.table <- function(otu, taxa, sort = TRUE){
  # ---- create otu table of lowest level
  
  # include label for "not assigned" reads to lowest level
  tnames <- rownames(taxa)
  taxa[which(taxa[,2] %in% "Not assigned"), length(taxa)] <- "Not_assigned "
  taxa <- data.frame(taxa[,ncol(taxa)])
  rownames(taxa) <- tnames
  colnames(taxa) <- "description"
  
  otu.table <- data.frame(merge(taxa, otu, by='row.names'), stringsAsFactors = FALSE)
  
  # remove rows with counts from higher level
  otu.table <- otu.table[!is.na(otu.table$description),]
  # set row name to NCBI ID
  rownames(otu.table) <-  stringr::str_extract(otu.table$description, "^(\\w+)")
  # remove taxa and ID rw
  otu.table$Row.names <- NULL
  # edit and rename description of NCBI ID
  otu.table$description <- sub("^(\\w+)\\s?(.*)$","\\2",otu.table$description )
  
  if (sort){
    otu.table <- otu.table[order(rownames(otu.table)),] # sort 
  }
  
  return(otu.table)
}
  
# create table from MEGAN6 biom file export, similar to QUIIME OTU table 
table2QUIIME.table <- function(otu, taxa){
  
  ## ---- generate QUIIME OTU tabe output 
  Q.table <- merge( format(otu, nsmall=1, trim=T), data.frame(apply(taxa,1, function(x) paste(x,collapse = "; ")))                   
                    , by='row.names') # match taxa to counts, set decimal
  colnames(Q.table)[c(1,length(colnames(Q.table)))] <- c('#OTU ID', 'taxonomy') # rename ID and taxa column 
  
  return(Q.table)

}


# function to rank OTU table by counts and report the featuede names of position
# min.lim to max.lim 
count_ranked_range <- function(otu, min.lim = -1, max.lim = -1){
  # get OTU sorted by mean ranking 
  
  # rank every column, subtract minimum
  otu.rank <- apply(otu, 2, function(x) rank(x)-min(rank(x)))
  # get rows from subsection
  otu.rank <- rownames(otu.rank[which( rownames(otu.rank) %in% names(sort(rowSums(otu.rank),
                                                                          decreasing=T))[min.lim:max.lim]),])
  return(otu.rank)
  
}


# function to exclude features from OTU while keeping the original sample size
exclude_rare <- function(x, names.stay){
  
  for (i in names(x)){
    # store column sum 
    otu.sub.sum <- colSums(x[[i]])
    # subset selection: components 
    x[[i]] <- x[[i]][names.stay, ]
    # restore original size 
    if ("Not_assigned" %in% names.stay){
      # join removed counts to not assigend 
      x[[i]]["Not_assigned",] <- x[[i]]["Not_assigned",] + otu.sub.sum - colSums(x[[i]])
    }else(
      # add new feature to OTU.sub with removed counts
      x[[i]][nrow(x[[i]])+1,] <- list("Not_assigned" = otu.sub.sum - colSums(x[[i]]))
    )
  }
  return(x)
}

exclude_rare_QUIIME <- function(otu.q, names.stay = otu.rank, resetID = T){
  
  # get OTU Lable from taxonomy
  taxon <- data.frame("#OTU ID"= stringr::str_split_fixed(
  stringr::str_split_fixed(otu.q$taxonomy, "; ", 4)[,4], " ", 2)[,1], 
  check.names = F, stringsAsFactors = F)
  rownames(taxon) <- rownames(otu.q)
  # get location of "Not_assigned"
  na.row <- which(rownames(taxon)== "-2")
  # get columns with counts vectro
  c.otu <-1:(ncol(otu.q)-1)
  # set na-label in taxon (for compare)
  taxon[na.row,1] <- "Not_assigned"
  # store sum per column
  otu.q.colsum <- colSums(otu.q[,c.otu])
  # set all rare species columns to zero 
  otu.q[which(!(taxon[,1] %in% otu.rank)), c.otu ] <- 0
  # restore sample sizes
  otu.q[na.row, c.otu] <- otu.q[na.row, c.otu] + otu.q.colsum - colSums( otu.q[ ,c.otu])
  
  
  # reset row names
  if (resetID){
    taxon$`#OTU ID`[which(taxon$`#OTU ID`== "NA")] <- which(taxon$`#OTU ID`== "NA")
    rownames(otu.q) <- taxon$`#OTU ID`
  }
  return(otu.q)
  
}










  