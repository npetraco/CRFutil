#' Quick generate a lattice graph
#'
#' Quick generate a lattice graph
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.lattice <- function(num.rows, num.cols, node.names.vec=NULL, cross.linksQ=T, plotQ=F){

  if(num.rows <= 1){
    stop("num.rows > 1")
  }

  if(num.cols <= 1){
    stop("num.cols > 1")
  }

  if(is.null(node.names.vec)){
    node.names <- paste0("X.",1:(num.rows*num.cols), sep="")
  } else {
    node.names <- node.names.vec
  }
  #print(node.names)

  jmat <- t(array(node.names ,c(num.rows, num.cols)))
  jmat

  eq <- NULL

  # Horizontal links
  for(i in 1:nrow(jmat)){
    for(j in 1:(ncol(jmat)-1)){
      htrm <- paste0(jmat[i,j], ":", jmat[i,j+1])
      #print(htrm)
      eq <- c(eq, htrm)
    }
  }

  if(cross.linksQ == T) {
    #Cross links and vertical links
    for(i in 1:(nrow(jmat)-1)) {
      for(j in 1:ncol(jmat)){
        if((j != 1)&(j != ncol(jmat))){
          #print(paste(jmat[i,j], "is NOT a corner"))
          #print(paste("   ", jmat[i+1,j-1], "is a neighbor"))
          #print(paste("   ", jmat[i+1,j],   "is a neighbor"))
          #print(paste("   ", jmat[i+1,j+1], "is a neighbor"))

          cterm1 <- paste0(jmat[i,j], ":", jmat[i+1,j-1])
          cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
          cterm3 <- paste0(jmat[i,j], ":", jmat[i+1,j+1])
          #print(cterm1)
          #print(cterm2)
          #print(cterm3)
          eq <- c(eq, cterm1, cterm2, cterm3)
        } else {
          #print(paste(jmat[i,j], "is a corner"))
          if(j == 1) {
            #print("    Start corner")
            #print(paste("   ", jmat[i+1,j],   "is a neighbor"))
            #print(paste("   ", jmat[i+1,j+1], "is a neighbor"))

            cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
            cterm3 <- paste0(jmat[i,j], ":", jmat[i+1,j+1])
            #print(cterm2)
            #print(cterm3)
            eq <- c(eq, cterm2, cterm3)
          } else {
            #print("    Stop corner")
            #print(paste("   ", jmat[i+1,j-1], "is a neighbor"))
            #print(paste("   ", jmat[i+1,j],   "is a neighbor"))

            cterm1 <- paste0(jmat[i,j], ":", jmat[i+1,j-1])
            cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
            #print(cterm1)
            #print(cterm2)
            eq <- c(eq, cterm1, cterm2)
          }
        }
      }
      #print("---------")
    }

  } else {
    # Vertical links
    for(j in 1:ncol(jmat)){
      for(i in 1:(nrow(jmat)-1)) {
        vtrm <- paste0(jmat[i,j], ":", jmat[i+1,j])
        #print(vtrm)
        eq <- c(eq, vtrm)
      }
    }

  }

  #print(eq)
  the.grphf.eq <- paste0("~", eq[1], sep="")
  for(i in 2:length(eq)){
    the.grphf.eq <- paste0(the.grphf.eq, "+", eq[i], sep="")
  }
  the.grphf.eq <- as.formula(the.grphf.eq)

  # I DONT KNOW WHY THIS ISNT WORKING
  # if(plotQ==TRUE){
  #   if(!is.null(dev.list())){
  #     dev.off()
  #     print("Here!")
  #   }
  #
  #   gpp <- ug(the.grphf.eq, result = "graph")
  #   plot(gpp)
  #   #plot(ug(the.grphf.eq, result = "graph"))
  # }

  return(the.grphf.eq)

}
