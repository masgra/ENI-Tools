


write_th_file = function(x, methods, path.out){
  # Iterate reactors
  for (j in 1:length(x[[1]])) {
    # set reactor th-file name
    fileConn <- file(paste(path.out,names(x[[1]])[j],'.txt',sep=''))
    th.lines <- c()
    # write th values to file
    for (i in 1:length(x)){
      th.lines <- c(th.lines, 
                    paste(methods[i],'~upperThreshold=',round(unlist(x[[i]][j])[1],digits=4),sep=''),
                    paste(methods[i],'~lowerThreshold=',round(unlist(x[[i]][j])[2],digits=4),sep='') )
    }
    # wirite emty line
    writeLines(th.lines,fileConn)
    # close fiel
    close(fileConn)
  }
}


weight_filtering = function(data , p, how="rank1"){
  # selects top and bottom percentil/ranked entries of a matrix and sets each other value of the matrix to zero.
  # "rank1" ranks absolute values
  # "rank2" ranks upper and lower values separately 
  # CAUTION: no NAN handling implemented yet
  
  if (how=="rank1"){
    # one sidet ranking
    data[!(rank(-abs(data), na.last=T, ties.method = "random")< p) ] <- 0
  }else if(how=="rank2"){
    # two sidet ranking
    data[!((rank(data, na.last=T, ties.method = "random")<round(p/2)) | (rank(-data, na.last=T, ties.method = "random")<floor(p/2)))] <- 0
  }else if (how == "percentil"){
    p <-.5 - abs(p-0.5) 
    lim <- quantile(data, c(percentil, 1-p))
    data[which((data>lim[1] & data<lim[2]))] <- 0
    
  }else{
    stop("'how' was not specified correctly. Options are 'rank1' (default), 'rank2' and 'percentil'.")
  }
  return(data)
}



