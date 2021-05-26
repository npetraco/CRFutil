#' Log sum exp trick. From Brendon Brewer's DNest code:
#'
#' Handy for calculating Z from a vector of log potentials.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
logsumexp <- function(logv) {
  n <- length(logv)
  max.logv <- max(logv)

  answer <- 0

  for(i in 1:n){
    answer <- answer + exp(logv[i] - max.logv)
  }
  answer <- max.logv + log(answer);

  return(answer)

}


#' Log sum exp trick. From Brendon Brewer's DNest code:
#'
#' Handy for calculating Z from a vector of log potentials.
#' A less readable but shorter Log sum exp:
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
logsumexp2 <- function(logv)
{
  n <- length(logv)
  max.logv <- max(logv)

  answer <-  max.logv + log(cumsum(c(0,exp(logv - max.logv)))[n+1])

  return(answer)

}


#' Code from prodlim library to match a row in a matrix
#'
#' XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
row.match <- function (x, table, nomatch = NA)   # **********NEEDS TO BE C
{
  if (class(table) == "matrix")
    table <- as.data.frame(table)
  if (is.null(dim(x)))
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}

#' Spits out permutation to re-order configs in targ with respect to ref
#'
#' XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
reorder_configs <- function (ref,targ){  # NEEDS TO BE C

  # FIX: When config is found take it out so next search is shorter.
  reord.idxs <- sapply(1:nrow(ref), function(xx){row.match(x = ref[xx,], table = targ)})

  # reordr.idxs       <- array(NA, nrow(ref))
  # targ.running      <- targ                 # This will shrink as reference configs are found
  # targ.running.idxs <- 1:nrow(targ)         # To keep track of which congis found and removed
  #
  # for(i in 1:nrow(ref)){
  #
  #   print(i)
  #
  #   # Where is the ref in the current target table?:
  #   tr.idx         <- row.match(x = ref[i,], table = targ.running)
  #
  #   # Translate this index to what it was in the original target table
  #   reordr.idxs[i] <- targ.running.idxs[tr.idx]
  #
  #   # Take out of the running what was found:
  #   targ.running      <- targ.running[-tr.idx, ]
  #   targ.running.idxs <- targ.running.idxs[-tr.idx]
  #
  # }

  return(reord.idxs)

}


#' Convenience function to complement a node in a configuration.
#' Assumes CRF 1,2 states
#'
#' XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
complement.at.idx <- function(configuration, complement.index){

  new.configuration <- configuration
  if(new.configuration[complement.index] == 1) {
    new.configuration[complement.index] <- 2
  } else if(new.configuration[complement.index] == 2){
    new.configuration[complement.index] <- 1
  } else {
    stop("States must be CRF states, i.e. 1 or 2!")
  }

  return(new.configuration)

}


#' Function to find which nodes are associated with which parameters.
#' XXXX
#'
#' The list produced maps node indices to the parameters they are asociated with
#'
#' The function will XXXX
#'
#' @param crf The XX
#' @param storeQ Logical, whether or not to store the the node to parameter list in the crf object
#' @return The function will XX
#'
#'
#' @export
nodes2params.list <- function(crf, storeQ = FALSE){

  node.par.assoc <- rep(list(NULL), crf$n.nodes)

  # Find the parameters associated with each node:
  for(i in 1:crf$n.nodes) {
    par.idxs  <- as.numeric(crf$node.par[i,,]) # parameters for the node

    # Get rid of the 0s:
    zeros.idxs <- which(par.idxs == 0)
    if(length(zeros.idxs) != 0) {
      par.idxs  <- par.idxs[-zeros.idxs]
    }

    if(length(par.idxs) == 0){
      stop("No parameters found for node: ", i) # **** This can happen when nodes or edges are not assigned parameters.
    }

    node.par.assoc[[i]] <- c(node.par.assoc[[i]], par.idxs)
    node.par.assoc[[i]] <- unique(node.par.assoc[[i]])
  }

  # Find the parameters associated with the nodes of each edge:
  for(i in 1:length(crf$edge.par)) {
    node.idx1 <- crf$edges[i,1] # node 1 of edge
    node.idx2 <- crf$edges[i,2] # node 2 of edge
    par.idxs  <- as.numeric(crf$edge.par[[i]]) # parameters for the edge

    # Get rid of the 0s:
    zeros.idxs <- which(par.idxs == 0)
    if(length(zeros.idxs) != 0) {
      par.idxs  <- par.idxs[-zeros.idxs]
    }

    if(length(par.idxs) == 0){
      stop("No parameters found for edge: ", i) # **** Is this ever allowed????
    }

    node.par.assoc[[node.idx1]] <- c(node.par.assoc[[node.idx1]], par.idxs)
    node.par.assoc[[node.idx2]] <- c(node.par.assoc[[node.idx2]], par.idxs)
    node.par.assoc[[node.idx1]] <- unique(node.par.assoc[[node.idx1]])
    node.par.assoc[[node.idx2]] <- unique(node.par.assoc[[node.idx2]])
  }

  if(storeQ == TRUE) {
    crf$nodes2pars <- node.par.assoc
  }

  return(node.par.assoc)
}


#' Function to find which nodes are associated with which parameters.
#' Generalization experiment
#'
#' The list produced maps node indices to the parameters they are asociated with
#'
#' The function will XXXX
#'
#' @param crf The XX
#' @param storeQ Logical, whether or not to store the the node to parameter list in the crf object
#' @return The function will XX
#'
#'
#' @export
nodes2params.list2 <- function(crf, storeQ = FALSE){

  node.par.assoc <- rep(list(NULL), crf$n.nodes)

  # Find the parameters associated with each node:
  for(i in 1:crf$n.nodes) {
    par.idxs  <- as.numeric(crf$node.par[i,,]) # parameters for the node

    # Get rid of the 0s:
    zeros.idxs <- which(par.idxs == 0)
    if(length(zeros.idxs) != 0) {
      par.idxs  <- par.idxs[-zeros.idxs]
    }

    if(length(par.idxs) == 0){
      warning("No parameters found for node: ", i) # **** This can happen when nodes or edges are not assigned parameters.
    } else {
      node.par.assoc[[i]] <- c(node.par.assoc[[i]], par.idxs)
      node.par.assoc[[i]] <- unique(node.par.assoc[[i]])
    }

  }

  # Find the parameters associated with the nodes of each edge:
  for(i in 1:length(crf$edge.par)) {
    node.idx1 <- crf$edges[i,1] # node 1 of edge
    node.idx2 <- crf$edges[i,2] # node 2 of edge
    par.idxs  <- as.numeric(crf$edge.par[[i]]) # parameters for the edge

    # Get rid of the 0s:
    zeros.idxs <- which(par.idxs == 0)
    if(length(zeros.idxs) != 0) {
      par.idxs  <- par.idxs[-zeros.idxs]
    }

    if(length(par.idxs) == 0){
      warning("No parameters found for edge: ", i) # **** Is this ever allowed????
    } else {
      node.par.assoc[[node.idx1]] <- c(node.par.assoc[[node.idx1]], par.idxs)
      node.par.assoc[[node.idx2]] <- c(node.par.assoc[[node.idx2]], par.idxs)
      node.par.assoc[[node.idx1]] <- unique(node.par.assoc[[node.idx1]])
      node.par.assoc[[node.idx2]] <- unique(node.par.assoc[[node.idx2]])
    }

  }

  if(storeQ == TRUE) {
    crf$nodes2pars <- node.par.assoc
  }

  return(node.par.assoc)
}


#' Function to find which parameters are associated with which nodes.
#' XXXX
#'
#' The list produced maps parameters indices to the nodes they are asociated with
#'
#' The function will XXXX
#'
#' @param crf The XX
#' @param storeQ Logical, whether or not to store the the parameter to node list in the crf object
#' @return The function will XX
#'
#'
#' @export
params2nodes.list <- function(crf, storeQ = FALSE){

  if(is.null(crf$nodes2par)){
    nodes2par <- nodes2params.list(crf, storeQ = FALSE)
  } else {
    nodes2par <- crf$nodes2par
  }

  par.node.assoc <- rep(list(NULL), crf$n.par)
  for(i in 1:crf$n.nodes) {
    for(k in 1:length(nodes2par[[i]])) {
      par.idx <- nodes2par[[i]][k]
      par.node.assoc[[par.idx]] <- c(par.node.assoc[[par.idx]] , i)
      par.node.assoc[[par.idx]] <- unique(par.node.assoc[[par.idx]])
    }
  }

  if(storeQ == TRUE) {
    crf$pars2nodes <- par.node.assoc
  }

  return(par.node.assoc)

}


#' Form adjacency matrix from edge matrix.
#' XXXX
#'
#' Assumes edges/nodes are numeric.
#'
#' The function will XXXX
#'
#' @param edge.mat The XX
#' @return The function will XX
#'
#'
#' @export
edges2adj <- function(edge.mat, plotQ=FALSE){

  num.nods <- max(edge.mat)
  adj.mat <- array(0, c(num.nods,num.nods))

  for(i in 1:nrow(edge.mat)){
    adj.mat[edge.mat[i,1], edge.mat[i,2]] <- 1
    adj.mat[edge.mat[i,2], edge.mat[i,1]] <- 1
  }
  colnames(adj.mat) <- 1:num.nods # **** Assumes edges/nodes are numeric and in order.
  rownames(adj.mat) <- 1:num.nods

  if(plotQ==TRUE){
    new.gph <- as(adj.mat,"graphNEL")
    if(!is.null(dev.list())){
      dev.off()
    }
    iplot(new.gph)
  }

  return(adj.mat)

}

#' Form pair-wise formula from adjacency matrix.
#' XXXX
#'
#' Assumes edges/nodes are numeric. Uses row and col numbers as node names.
#'
#' The function will XXXX
#'
#' @param edge.mat The XX
#' @return The function will XX
#'
#'
#' @export
adj2formula <- function(adj.mat, Xoption=FALSE){

  # To be safe, just use row and col numbers as node names for now. Add X if required by a
  # formula interface

  num.nodes <- ncol(adj.mat)
  if(Xoption == TRUE) {
    node.idxs <- apply(as.matrix(paste0("X",1:num.nodes)), 2, paste0, collapse = " + ")
  } else {
    node.idxs <- apply(as.matrix(1:num.nodes), 2, paste0, collapse = " + ")
  }

  #print(node.idxs)

  adj.up.tri <- adj.mat*upper.tri(adj.mat)
  edgs.idxs  <- which(adj.up.tri==1, arr.ind = T)

  # Order edges by first node
  edg.ord   <- order(edgs.idxs[,1])
  #edgs.idxs <- edgs.idxs[edg.ord,]

  if(Xoption == TRUE){
    edgs.idxs <- cbind(paste0("X",edgs.idxs[edg.ord,1]),paste0("X",edgs.idxs[edg.ord,2]))
  } else {
    edgs.idxs <- edgs.idxs[edg.ord,]
  }
  #print(edgs.idxs)

  edgs.idxs <- apply(edgs.idxs, 1, paste0, collapse = ":")
  edgs.idxs <- apply(as.matrix(edgs.idxs), 2, paste0, collapse = " + ")
  #print(edgs.idxs)

  graph.fomla <- as.formula(paste0("~", node.idxs," + ", edgs.idxs))

  return(graph.fomla)

}

