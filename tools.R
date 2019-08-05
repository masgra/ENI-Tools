ReadCountTable <- function(files , header=TRUE, sep =" " , dec = "." , row.names = 1, col.names= TRUE,
                            check.names = FALSE, skip = 0, nrows = -1, fill = TRUE,
                            re.ID.separator = NULL , ID.split.elem = 1, 
                            re.samples = NULL){
  ## reads in a list files with count data and joins dem by their row names
  # optional: simplifies the row names
  # optional: aggregates the columns by their names (by sum)
  
  # Args:
  #   files: list of data file paths
  #   header, sep, dec, row.names, col.names, check.names, skip, nrows, fill: see csv.read()
  #   re.ID.separator: separator to identify the row name ID.
  #   ID.split.elem: element to use as row name ID (usually the first)
  #   re.samples: regular expression for sample name
  #
  # Returns:
  #   ds: merged dataset of all records
  
  ds <- data.frame()
  
  # load and process each file and merge it with final table
  for (i in 1:length(files)){
    
    # load data from files
    ds.i <-read.csv(files[i], header = header, sep = sep, dec = dec,
                    row.names = row.names, check.names = check.names, 
                    skip = skip, nrows = nrows, fill = fill)
    
    # edit row names: only keep the regular expression
    if (!is.null(re.ID.separator)){
      rownames(ds.i) <- lapply(rownames(ds.i), function(x) unlist(strsplit(x,re.ID.separator))[ID.split.elem])
    }
    
    # merge all fields by row names
    ds <- merge(ds,ds.i, all=T, by="row.names")
    rownames(ds) <-ds$Row.names
    ds <- ds[,-1]
  }
  
  # aggregate records by column names (by the output of a regular expression )
  if (!is.null(re.samples)){
    samples <- sub(re.samples, "\\1", colnames(ds),perl=TRUE)
    ds <- aggregate( subset(data.frame(t(ds))), 
                     by=list(x.names = samples), FUN='sum' ) 
    rownames(ds) <- ds$x.names
    ds <-data.frame(t(subset(ds, select = -x.names)),check.names = FALSE)
  }
  # set NaN to zero 
  ds[is.na(ds)] <- 0
  return(ds)
}



Biome2Table <- function(filepath, sample.reg=NULL, rename=NULL, sort.rows=TRUE, sort.cols=TRUE ){ 
  ## reads in a biom file and joins thier records by matching row names
  # optional: simplifies the row names
  # optional: aggregates the columns by their names (by sum)
  
  # Args:
  #   filepath: file path of .biom file
  #   sample.reg: regular expression to identify duplicated records (otherwise no join)
  #   renam: regular expression of string that should be kept as sample ID 
  #   sort.rows: flag to sort rows alphabetically (strings) or numerically (numerics)
  #   sort.cols: flag to sort columns alphabetically (strings) or numerically (numerics)
  #   
  #
  # Returns:
  #   otu: (joined) dataset of all records, with ortholog
  
  
  # load data form biom-format file
  data <- phyloseq::import_biom(filepath)
  # get OTU table as data frame
  otu <- data.frame(data@otu_table@.Data, check.names=F, stringsAsFactors = FALSE)
  # get ortholog, levels separated as columns
  ortholog <-data.frame(data@tax_table@.Data, stringsAsFactors = FALSE)
  if(!is.null(sample.reg)){
    # sum samples with same ID
    otu <- sapply( unique( stringr::str_extract(colnames(otu), sample.reg) ), 
                   function(x) rowSums(otu[,grep(x, colnames(otu)), drop=FALSE]))
  }
  # rearrange order of columns to increasing by sample number  
  if(!is.null(rename)){
    colnames(otu) <- stringr::str_extract(colnames(otu), "[0-9]+") # rename
      
  }
  # sort by record names
  if (sort.cols){
    if (!any(is.na(as.numeric(colnames(otu))))){
      # sort numerically if record names are numbers
      otu <- otu[,order(as.numeric(colnames(otu)))] # sort
    }else{
      # otherwise sort records alphabetically
      otu <- otu[,order(colnames(otu))]
    }
  }
  
  otu <- cbind(otu,ortholog)
  # sort components 
  if (sort.rows){
    if (!any(is.na(as.numeric(rownames(otu))))){
      # sort numerically if component names are numbers
      otu <- otu[order(as.numeric(rownames(otu))),] # sort
    }else{
      # otherwise sort components alphabetically
      otu <- otu[order(rownames(otu)),]
    }
  }
  
  return (otu)
}








SelectRankedComponents <- function(otu, min.lim = -1, max.lim = -1){
  # This function computes the rank of each component per column, computes the average rank sum per component and 
  # returns the names of the components between the upper and lower average ranks.
  #
  # Args:
  #  otu: imput table with records as columns (rows=components get ranks)
  #  min.lim, max.lim: total average rank limits
  #
  # Return: 
  #  otu.rank: a name vector of the components within the defined rank range 
  
  # rank every column, subtract minimum (this reverses the order of ranks)
  otu.rank <- apply(otu, 2, function(x) rank(x)-min(rank(x)))
  # compute average sum of ranks per row, sort them and get row names of range subsection 
  otu.rank <- rownames(otu.rank[which( rownames(otu.rank) %in% names(sort(rowSums(otu.rank),
                                                                          decreasing=T))[min.lim:max.lim]),])
  return(otu.rank)
}







excludeComponents <- function(x0, names.stay, trash.node.name="Not_assigned"){
  # This function excludes all components from the count table that are not in the given "name.stay" vector.
  # Counts are moved to a trash node to preserve the closure. 
  #
  # Args:
  #  x0: imput table with records as columns (rows=components get ranks)
  #  names.stay: names vector with component names that will stay
  #  trash.node.name: name of trash node. A new one will be created if this row not already exists
  #
  # Return: 
  #  otu.rank: a name vector of the components within the defined rank range 
  
  
  
  # store original count table
  x <- x0[names.stay,]
  
  # create trashnode if not exists 
  if(sum(rownames(x) %in% trash.node.name) <= 0){
    # add trashnode if it is not part of the selected range
    x[nrow(x)+1,] <- 0
    rownames(x)[nrow(x)] <- trash.node.name
  }
  
  # get location of trash node
  trash.row <- which(rownames(x)== trash.node.name)
  
  # store exclude counts in trash node 
  x[trash.row, ] <- x[trash.row, ] + colSums(x0, na.rm = T) - colSums(x, na.rm = T)
  return(x)
  
}





writeThresholdFile = function(x, methods =NULL, path.out){
  # This function takes a list of threshold object as its input and writes them into a file format,
  # that can be read in form CoNet as an input argument. Multiple files will be written If the list 
  # holds thresholds for multiple experiments. 
  #
  # Args:
  #  x: input threshold list. Should look like: 
  #     list( metric1=list( experiment1=[upper bound, lower bound], 
  #                         experiment2 = [...],
  #                         ... ), 
  #           metric2 = ... ), ...)
  #  methods: optional - if names in the threshold list are different from the CoNet metric arguments
  #  path.out: place to store the output file(s) 
  #
  # Return: 
  #  none 
  
  # Iterate over experiments
  for (j in 1:length(x[[1]])) {
    # set experiment th-file name
    fileConn <- file(paste(path.out,names(x[[1]])[j],'.txt',sep=''))
    th.lines <- c()
    # write threshold values to file
    for (i in 1:length(x)){
      if (is.null(methods)){
        th.lines <- c(th.lines, 
                    paste(names(x)[i],'~upperThreshold=',round(unlist(x[[i]][j])[1],digits=4),sep=''),
                    paste(names(x)[i],'~lowerThreshold=',round(unlist(x[[i]][j])[2],digits=4),sep='') )
      } else if (length(methods)==length(x)){
        th.lines <- c(th.lines, 
                      paste(methods[i],'~upperThreshold=',round(unlist(x[[i]][j])[1],digits=4),sep=''),
                      paste(methods[i],'~lowerThreshold=',round(unlist(x[[i]][j])[2],digits=4),sep='') )
      }else{
        stop("size of input list (= ", length(x), ") in not equal to lenght of method vector (= ", length(methods), ")." )
      }
    }
    # write empty line
    writeLines(th.lines,fileConn)
    # close file
    close(fileConn)
  }
}






weightFiltering = function(data , p, how="rank1"){
  # selects top and bottom percentile/ranked entries of a matrix and sets each other value of the matrix to zero.
  # "rank1" ranks absolute values
  # "rank2" ranks upper and lower values separately 
  # CAUTION: no NAN handling implemented yet
  
  if (how=="rank1"){
    # one sidet ranking
    data[!(rank(-abs(data), na.last=T, ties.method = "random")<= 2*p) ] <- 0
  }else if(how=="rank2"){
    # two sidet ranking
    data[!((rank(data, na.last=T, ties.method = "random")<=p) | (rank(-data, na.last=T, ties.method = "random")<=p))] <- 0
  }else if (how == "percentile"){
    p <-.5 - abs(p-0.5) 
    lim <- quantile(data, c(p, 1-p))
    data[which((data>lim[1] & data<lim[2]))] <- 0
    
  }else{
    stop("'how' was not specified correctly. Options are 'rank1' (default), 'rank2' and 'percentile'.")
  }
  return(data)
}





getCytoscapeEdgeDS <- function(x, names, expand.Unreported.Verteces=TRUE, feature = NULL){
  # creates a data structure from a adjacency matrix input that can easily exported as an edge table for Cytoscape. 
  # unreported notes can be added at the end of the table. 
  # features for nodes can be added as well. 
  
  #
  # Args:
  #  x: adjacency matrix
  #  names: names of components (because SpiecEasi does not provide any names)
  #  expand.Unreported.Verteces: flag, whether to add components with no edges to the end of the structure
  #  feature: optional additional feature list. Names must match the names of components
  #
  # Return: 
  #  res: data frame of edges with weights, interactions and associations
  #      component1 (V1), Component2 (V2), weight, interaction, features ...  
  

  
  
  # set upper triangula to zero (to remove reverse interaction), set row and col name 
  w <- data.frame(x$th.weight*lower.tri(x$th.weight,diag=F), row.names = names)
  
  #colnames(w) <- names
  colnames(w) <- names
  # create melted table with each row as a possible interaction  
  res <- data.frame( reshape::melt(as.matrix(w)) )
  
  
  colnames(res) = c("V1", "V2", "weight")
  # remove rows with zero weights
  res <- res[res$weight!=0, ] 
  
  # create interaction item
  res$interaction <- getInteraction( x= c(res$weight) )
  
  if (expand.Unreported.Verteces){
    # add vertex name with no edges to ensure that all vertices’ are later plotted in cytoscape 
    r <- data.frame(names[which(!(names %in% as.character(res$V1)))])
    if(nrow(r)>0){
      r[,2:4] <- NA
      colnames(r) = c("V1", "V2", "weight","interaction")
      res <-rbind(res,r)
    }
  }
  if(!is.null(feature)){
    res <- cbind(res, feature[match(res$V1,rownames(feature)),])
  }
  return(res)
}


