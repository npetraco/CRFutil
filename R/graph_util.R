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


#' Send in two or more models in the form of a list of crf objects or edge matrices and compare edges
#'
#' Send in two or more models in the form of a list of crf objects or edge matrices and compare edges
#'
#' Assumes all models have the same number of nodes and node names are just numbers
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
compare_edges <- function(model.list, num.nodes=NULL){

  num.models <- length(model.list)

  model.type <- class(model.list[[1]]) # Look at the first element of the list and determine the data's type
  if(model.type == "CRF") {
    print("Models are CRF objects.")

    crf.obj.loc           <- model.list[[1]]
    num.nodes.loc         <- crf.obj.loc$n.nodes
    saturated.edge.matrix <- array(-1, c(num.nodes.loc*(num.nodes.loc-1)/2, 2))
    print(dim(saturated.edge.matrix))

    # Enumerate all possible edges. This will be the reference to compare models
    count <- 1
    for(i in 1:num.nodes.loc) {
      for(j in 1:num.nodes.loc) {
        if(i < j) {
          #print(count)
          saturated.edge.matrix[count, 1] <- i
          saturated.edge.matrix[count, 2] <- j
          count <- count + 1
        }
      }
    }
    #print(saturated.edge.matrix)

    edgeQ.mat <- NULL
    for(i in 1:num.models){

      crf.obj.loc <- model.list[[i]]
      edge.mat.loc <- crf.obj.loc$edges
      #print(edge.mat.loc)

      #edgeQ.vec <- numeric(nrow(saturated.edge.matrix))
      edgeQ.vec <- array(-1, nrow(saturated.edge.matrix))
      for(j in 1:nrow(saturated.edge.matrix)) {
        #edgeQ.vec[j] <- as.numeric( ((saturated.edge.matrix[j,1] %in% edge.mat.loc[,1]) & (saturated.edge.matrix[j,2] %in% edge.mat.loc[,2]) ) | (saturated.edge.matrix[j,2] %in% edge.mat.loc[,1]) & (saturated.edge.matrix[j,1] %in% edge.mat.loc[,2]) )
        # Number of times (possible) edge observed in edge matrix of graph. Should be 0 or 1.
        # If two or more, throw an error.
        num.times.edge.obs <- sum(sapply(1:nrow(edge.mat.loc), function(xx){sum(sort(saturated.edge.matrix[j,]) == sort(edge.mat.loc[xx,])) == 2}))
        if(num.times.edge.obs >= 2) {
          print(saturated.edge.matrix[j,])
          stop(paste0("Edge above appears in edge matrix of graph ", i, "  more than once. Something is wrong!"))
        }

        edgeQ.vec[j] <- num.times.edge.obs # Should be 0 or 1 at this point

      }
      #print(cbind(saturated.edge.matrix, edgeQ.vec))
      edgeQ.mat <- cbind(edgeQ.mat, edgeQ.vec)

    }

  } else if( ("matrix" %in% model.type) | ("array" %in% model.type) ) {
    print("Models are edge matrices.")

    # *********Need num.nodes for these

  } else {
    stop("model.list must be a list of CRF objects or edge matrices!")
  }

  print(data.frame(saturated.edge.matrix, edgeQ.mat))

}
