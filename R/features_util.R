#' Utility function to compute unconditional features phi_i(config)
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
phi.features <- function(config, edges.mat, node.par, edge.par, ff) {

  num.nodes  <- length(config)
  num.edges  <- nrow(edges.mat)
  num.params <- max(unlist(list(node.par, edge.par)))
  phi.vec    <- numeric(num.params)

  phi.idx <- 0

  # Nodes: \phi_i({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \tau}_i, 0}
  for(i in 1:num.nodes){
    phi.idx <- as.numeric(ff(config[i]) %*% node.par[i,,1])
    #print(paste("phi index for node ", i, " ", phi.idx))
    if(phi.idx != 0) {
      phi.vec[phi.idx] <- 1 # Don't do this???? Assumes parameters have unique consecutive indices... Causing SUBTLE BUG?????
      #phi.vec[i] <- 1 #???????????
    }
    #print("phi.vec")
    #print(phi.vec)
  }

  # Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  for(i in 1:num.edges) {

    # Don't do this?? Assumes parameters have unique consecutive indices...
    # # Check and see if we reached the end of phi. No point in doing the rest of the edges if we did:
    if(phi.idx == num.params) {
      break()
    }

    phi.idx <- as.numeric(ff(config[edges.mat[i,1]]) %*% edge.par[[i]][,,1] %*% ff(config[edges.mat[i,2]]))
    #print(paste("phi index for edge ", i, " ", phi.idx))
    if(phi.idx != 0) {
      phi.vec[phi.idx] <- 1
    }
  }

  return(phi.vec)
}


#' Utility function to compute labeled explicit unconditional features phi_i(config)
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
phi.features.explicit <- function(config, edges.mat, node.par, edge.par, ff) {

  num.nodes   <- length(config)
  num.edges   <- nrow(edges.mat)
  num.params  <- max(unlist(list(node.par, edge.par)))
  phi.vec     <- numeric(num.params)
  phi.vec.lbl <- numeric(num.params)

  phi.idx <- 0

  # Nodes: \phi_i({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \tau}_i, 0}
  for(i in 1:num.nodes){
    phi.idx              <- Eone(config[i], node.par[i,,1], ff) # Just re-use the Eone code
    phi.i                <- 1 - (phi.idx == 0)
    phi.vec[phi.idx]     <- phi.i
    phi.vec.lbl[phi.idx] <- phi.idx
  }

  # Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  for(i in 1:num.edges) {
    phi.idx              <- Etwo(config[edges.mat[i,1]], config[edges.mat[i,2]], edge.par[[i]][,,1], ff) # Just re-use the Etwo code
    phi.i                <- 1 - (phi.idx == 0)
    phi.vec[phi.idx]     <- phi.i
    phi.vec.lbl[phi.idx] <- phi.idx
  }

  return(list(phi.vec, phi.vec.lbl))
}

#' Utility function to convert node and edge indices to the
#' parameter index they are associated with.
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.par.idx <- function(config, i=NULL, j=NULL, node.par=NULL, edge.par=NULL, edge.mat=NULL, ff, printQ=FALSE) {

  if(!is.null(i) & !is.null(j)) {
    if(!is.null(node.par)) {
      stop("Two node indices specified but node.par entered. Enter ONLY edge.par with two node indices.")
    }
  }

  if(is.null(i) & is.null(j)){
    stop("No node indices entered!")
  } else if(!is.null(i)){
    if(!is.null(j)) { # Compute an edge parameter index
      if(is.null(edge.par)){
        stop("Two node indices entered but no edge.par input!")
      } else {
        if(is.null(edge.mat)) {
          stop("edge.par specified but there is no edge.mat!")
        } else {
          edge.idx <- row.match(c(i,j), edge.mat)
          if(length(edge.idx)==0){
            stop("Edge not found! Check i, j, and edge.mat entered.")
          } else {

            if(printQ==T){
              print(paste("Node i:", i))
              print(paste("Node j:", j))

              print("f0(i):")
              print(ff(config[i]))
              print(paste("class f0i:",class(ff(config[i]))))

              print("f0(j):")
              print(ff(config[j]))
              print(paste("class f0j:",class(ff(config[i]))))

              print(paste("Edge index: ",edge.idx))
              print("Edge:")
              print(edge.mat[edge.idx,])

              print("edge.par matrix for edge:")
              print(edge.par[[edge.idx]][,,1])
              print(paste("class: ", class(edge.par[[edge.idx]][,,1])))

              print("===================================")
            }

            par.idx <- as.numeric(ff(config[i]) %*% edge.par[[edge.idx]][,,1] %*% ff(config[j]))
          }
        }
      }
    } else { # Compute a node parameter index
      if(is.null(node.par)) {
        stop("One node indices entered but no node.par input!")
      } else {
        if(!is.null(edge.mat)){
          warning("Node parameter index calculation requested but an edge.mat is specified. Ignoring edge.mat.")
        }
        par.idx <- as.numeric(ff(config[i]) %*% node.par[i,,1])
      }
    }
  } else {
    stop("No node i index entered!")
  }

  return(par.idx)

}


#' Utility function to convert a parameter index to the node/edge indices it is associated with.
#'
#' Sort of an inverse function of get.par.idx. Intended primarily for checking/debugging.
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.node.idxs <- function(par.idx, node.par=NULL, edge.par=NULL, edge.mat=NULL) {
  node.idxs <- NULL
  edge.idxs <- NULL
  all.idxs  <- NULL

  if(is.null(edge.par)){
    warning("No edge.par specified.")
  }

  node.idxs <- c(which(node.par[,1,1]==par.idx), which(node.par[,2,1]==par.idx))
  for(i in 1:nrow(edge.mat)){
    #print(edge.par[[i]][,,1])
    if(par.idx %in% edge.par[[i]][,,1]) {
      #print("HERE")
      edge.idxs <- rbind(edge.idxs, edge.mat[i,])
    }
  }

  all.idxs <- list(node.idxs, edge.idxs) # Just unlist and unique to put into a vector
  names(all.idxs) <- c("node.idxs","edge.idxs")

  return(all.idxs)

}


#' Utility function to compute phi component CONDITIONED for specific node or edge
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function can be used to compute the CONDITINAL phi feature vector used for pseudolikelihood calculations.
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
phi.component <- function(config, i=NULL, j=NULL, node.par=NULL, edge.par=NULL, edge.mat=NULL, ff) {

  par.idx <- get.par.idx(config, i, j, node.par, edge.par, edge.mat, ff)
  comp <- 1 - (par.idx == 0)

  return(comp)

}


#' Utility function to compute model matrix
#'
#' Model matrix is the matrix of unconditional phi "features". Column
#' sums of the model matrix should be the sufficient statistics.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
compute.model.matrix <- function(configs, edges.mat, node.par, edge.par, ff) {

  model.mat <- t(sapply(1:nrow(configs), function(xx){
    phi.features(config = configs[xx,], edges.mat, node.par, edge.par, ff)
  }))

  return(model.mat)

}


#' Utility function to compute mean vector of unconditional features with supplied theta in a CRF network object.
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
feature.means <- function(crf, inference.func = infer.exact) {

  inference.info <- inference.func(crf)

  # Compute \text{E}_{\hat{\theta}_i}[{\phi_i}] with marginals at the optimal theta:
  # I.E. use "inference" to avoid computing X and Z directly. Use "junction tree":
  nodeMap    <- as.numeric(crf$node.par[,,]) # Flatten node parameter index matrix into a vector
  edgeMap    <- unlist(crf$edge.par)         # Flatten edge parameter index matries into a vector
  num.params <- max(c(nodeMap,edgeMap))      # eventually change to crf$n.par

  nodeBel    <- as.numeric(inference.info$node.bel) # Flatten node marginals matrix into a vector
  edgeBel    <- unlist(inference.info$edge.bel)     # Flatten edge marginals matrices into a vector

  param.phi.means <- numeric(num.params)
  for(i in 1:num.params) {
    if(i %in% nodeMap) {
      #print(paste("Param", i,"is a node parameter"))
      param.idxs <- which(nodeMap == i) # indices for parameter i
      param.phi.means[i] <- sum(nodeBel[param.idxs])
    }
    if(i %in% edgeMap) {
      #print(paste("Param", i,"is an edge parameter"))
      param.idxs <- which(edgeMap == i) # indices for parameter i
      param.phi.means[i] <- sum(edgeBel[param.idxs])
    }
  }

  return(param.phi.means)
}
