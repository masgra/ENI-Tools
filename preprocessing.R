


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
  

table2QUIIME.table <- function(otu, taxa){
  
  ## ---- generate QUIIME OTU tabe output 
  Q.table <- merge( format(otu, nsmall=1, trim=T), data.frame(apply(taxa,1, function(x) paste(x,collapse = "; ")))                   
                    , by='row.names') # match taxa to counts, set decimal
  colnames(Q.table)[c(1,length(colnames(Q.table)))] <- c('#OTU ID', 'taxonomy') # rename ID and taxa column 
  
  return(Q.table)

}




  