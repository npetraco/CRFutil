#' Utility function to compute features phi_i(config)
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
    if(phi.idx != 0) {
      phi.vec[phi.idx] <- 1
    }
  }

  # Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  for(i in 1:num.edges) {

    # Check and see if we reached the end of phi. No point in doing the rest of the edges if we did:
    if(phi.idx == num.params) {
      break()
    }

    phi.idx <- as.numeric(ff(config[edges.mat[i,1]]) %*% edge.par[[i]][,,1] %*% ff(config[edges.mat[i,2]]))
    if(phi.idx != 0) {
      phi.vec[phi.idx] <- 1
    }
  }

  return(phi.vec)
}


#' Utility function to compute labeled explicit features phi_i(config)
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
    #phi.idx              <- as.numeric(ff(config[i]) %*% node.par[i,,1])
    phi.idx              <- Eone(config[i], node.par[i,,1], ff) # Just re-use the Eone code
    phi.i                <- 1 - (phi.idx == 0)
    phi.vec[phi.idx]     <- phi.i
    phi.vec.lbl[phi.idx] <- phi.idx
    # if(phi.idx != 0) {
    #   phi.vec[phi.idx] <- 1
    # }
  }

  # Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  for(i in 1:num.edges) {

    # Check and see if we reached the end of phi. No point in doing the rest of the edges if we did:
    # if(phi.idx == num.params) {
    #   break()
    # }

    #phi.idx              <- as.numeric(ff(config[edges.mat[i,1]]) %*% edge.par[[i]][,,1] %*% ff(config[edges.mat[i,2]]))
    phi.idx              <- Etwo(config[edges.mat[i,1]], config[edges.mat[i,2]], edge.par[[i]][,,1], ff) # Just re-use the Etwo code
    phi.i                <- 1 - (phi.idx == 0)
    phi.vec[phi.idx]     <- phi.i
    phi.vec.lbl[phi.idx] <- phi.idx
    # if(phi.idx != 0) {
    #   phi.vec[phi.idx] <- 1
    # }
  }

  return(list(phi.vec, phi.vec.lbl))
}


#' Utility function to compute model matrix
#'
#' Model matrix is the matrix of phi "features". Column
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

  model.mat <- t(sapply(1:nrow(configs), function(xx){phi.features(config = configs[xx,], edges.mat, node.par, edge.par, ff)}))
  return(model.mat)

}


#' Utility function to compute mean vector of features with supplied theta in a CRF network object.
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
  # node.param.phi.means <- rowSums(inference.info$node.bel * (crf$node.par>0)[,,])
  nodeMap <- crf$node.par[,,]    # IS THERE A BETTER WAY TO DO THIS THAN NESTED LOOPING??
  # for(i in 1:nrow(nodeMap)) {
  #   for(j in 1:ncol(nodeMap)){
  #
  #   }
  # }

  # edge.param.phi.means <- numeric(crf$n.edges)
  # for(i in 1:crf$n.edges) {
  #   edge.param.phi.means[i] <- sum(inference.info$edge.bel[[i]] * (crf$edge.par[[i]]>0)[,,])
  # }
  # mean.phi.vec <- c(node.param.phi.means,edge.param.phi.means)
  #
  # return(mean.phi.vec)

}

