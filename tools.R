


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
    writeLines(c(th.lines,''),fileConn)
    # close fiel
    close(fileConn)
  }
}


