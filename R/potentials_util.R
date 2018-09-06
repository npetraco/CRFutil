#' Shift node potentials back to the originally fit parameter vector
#'
#' For testing purposes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
shift.pots <- function(crf) {

  parms <- crf$par

  # Loop over elements of parameter index matrix (also called nodeMap and edgeMap in UGM)
  # and put elements of exp(parms) where they belong

  node.par         <- crf$node.par[,,1]
  node.pot.shifted <- array(1, c(dim(node.par),1) )
  for(i in 1:nrow(node.par)) {
    for(j in 1:ncol(node.par)) {
      if(node.par[i,j] != 0) {
        par.idx <- node.par[i,j]
        node.pot.shifted[i,j,1] <- exp(parms[par.idx])
      }
    }
  }


  edge.par <- crf$edge.par
  edge.pot.shifted <- rep(list(array(1,c(2,2))), length(edge.par))
  for(k in 1:length(edge.par)) {
    for(i in 1:nrow(edge.par[[k]])) {
      for(j in 1:ncol(edge.par[[k]])) {
        if(edge.par[[k]][i,j,1] != 0) {
          par.idx <- edge.par[[k]][i,j,1]
          edge.pot.shifted[[k]][i,j] <- exp(parms[par.idx])
        }
      }
    }
  }

  crf$node.pot <- node.pot.shifted
  crf$edge.pot <- edge.pot.shifted
  print("Potentials regenerated from parameter vector via node.par and edge.par")

}


#' Make potentials from input parameter vector.
#'
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.pots <- function(parms, crf, rescaleQ=FALSE, replaceQ=FALSE, printQ=FALSE) {

  # Loop over elements of parameter index matrix (also called nodeMap and edgeMap in UGM)
  # and put elements of exp(parms) where they belong
  node.par         <- crf$node.par[,,1]
  node.pot.shifted <- array(1, c(dim(node.par),1) )
  for(i in 1:nrow(node.par)) {
    for(j in 1:ncol(node.par)) {
      if(node.par[i,j] != 0) {
        par.idx <- node.par[i,j]
        node.pot.shifted[i,j,1] <- exp(parms[par.idx])
      }
      if(rescaleQ==TRUE) {
        node.pot.shifted[i,,1] <- node.pot.shifted[i,,1]/max(node.pot.shifted[i,,1])
      }
    }
  }

  edge.par <- crf$edge.par
  edge.pot.shifted <- rep(list(array(1,c(2,2))), length(edge.par))
  for(k in 1:length(edge.par)) {
    for(i in 1:nrow(edge.par[[k]])) {
      for(j in 1:ncol(edge.par[[k]])) {
        if(edge.par[[k]][i,j,1] != 0) {
          par.idx <- edge.par[[k]][i,j,1]
          edge.pot.shifted[[k]][i,j] <- exp(parms[par.idx])
        }
      }
    }
    if(rescaleQ==TRUE){
      #print(paste(k, " ", max(edge.pot.shifted[[k]]) ))
      edge.pot.shifted[[k]] <- edge.pot.shifted[[k]]/max(edge.pot.shifted[[k]])
    }

  }

  if(replaceQ==TRUE){
    crf$node.pot <- node.pot.shifted
    crf$edge.pot <- edge.pot.shifted
    if(printQ==TRUE){
      print("Potentials computed from parameter vector, node.par and edge.par were written to CRF object's node.pot and edge.pot.")
    }
  }

  # Don't worry about the .shifted extension. It's just so I could copy code from shift.pots.
  return(list(node.pot.shifted,edge.pot.shifted))

}


#' Make a new parameter vector from ALL of the potential values contained in the crf object.
#' This makes new node/edge pars and returns them. If node/edge pars are contained in the crf object passed in,
#' they are ignored.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.par.from.all.potentials <- function(crf) {

  num.states  <- ncol(crf$node.pot)

  # Flatten node parameters into format:
  # tauA_s1, tauA_s2, tauA_s3, ... tauB_s1, tauB_s2, ....
  node.params <- as.numeric(t(crf$node.pot))
  node.par    <- t(array(seq(from=1, to=num.states*crf$n.nodes, by=1), c(ncol(crf$node.pot), nrow(crf$node.pot))))
  node.par    <- array(node.par, c(dim(node.par),1)) #Annoying, but CRF does this (think it was forgotten that it isn't necessary and never fixed...)

  # Flatten edge parameters into format:
  # omega_edge1_s1,s1, omega_edge1_s1,s2, omega_edge1_s1,s3, ... (matrix in list element 1 row stacked)
  # omega_edge1_s2,s1, omega_edge1_s2,s2, omega_edge1_s2,s3, ...
  # .
  # .
  # .
  # omega_edge2_s1,s1, omega_edge2_s1,s2, omega_edge2_s1,s3, ... (matrix in list element 2 row stacked)
  # .
  # .
  # .
  edge.par       <- rep(list(NULL), crf$n.edges)
  edge.params    <- NULL
  curr.param.num <- num.states*crf$n.nodes + 1 # Keep track of what param number we are on
  for(i in 1:crf$n.edges) {
    edge.i.params <- as.numeric(t(crf$edge.pot[[i]]))
    edge.params   <- c(edge.params, edge.i.params)

    edge.par[[i]] <- t(array(seq(from = curr.param.num, to = curr.param.num + length(edge.i.params), by = 1), dim(crf$edge.pot[[i]])))
    edge.par[[i]] <- array(edge.par[[i]], c(dim(edge.par[[i]]), 1))
    curr.param.num <- curr.param.num + length(edge.i.params)
  }

  new.params.vec <- log(c(node.params, edge.params))
  new.potentials.info <- list(new.params.vec, node.par, edge.par)
  names(new.potentials.info) <- c("par", "node.par", "edge.par")

  return(new.potentials.info)
}


#' This version makes a parameter vector from the potentials according to the node/edge pars contained in the
#' crf object. If the potentials are not scalled according to the codeing in node/edge par, they are re-scaled
#' before being used.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.par.from.potentials <- function(crf) {

  nodeMap <- crf$node.par
  edgeMap <- crf$edge.par
  nodePot <- crf$node.pot
  edgePot <- crf$edge.pot

  # First make sure that the potentials are scaled according to the node/edge pars:
  for(i in 1:nrow(nodePot)){

    zero.idxs <- which(nodeMap[i,,1] == 0)
    if(length(zero.idxs) == 1) {
      nodePot[i, ] <- nodePot[i,]/nodePot[i, zero.idxs] # These are pots, so divide not subtract
    }
    if(length(zero.idxs) > 2) {
      warning("More than 2 reference elements found in node.par row: ", i, ". Skipping that row.")
    }
    if(length(zero.idxs) == 0) {
      warning("No reference elements found in node.par row: ", i, ". Skipping that row.")
    }

  }
  #print(nodePot)

  for(k in 1:length(edgePot)) {
    for(i in 1:nrow(edgePot[[k]]) ) {

      zero.idxs <- which( edgeMap[[k]][i,,1] == 0 )
      if(length(zero.idxs) == 1) {
        edgePot[[k]][i, ] <- edgePot[[k]][i,]/edgePot[[k]][i, zero.idxs] # These are pots, so divide not subtract
      }
      if(length(zero.idxs) > 2) {
        warning("More than 2 reference elements found in edge.par: ", k, " row ", i, ". Skipping that row.")
      }
      if(length(zero.idxs) == 0) {
        warning("No reference elements found in edge.par: ", k, " row ", i, ". Skipping that row.")
      }
    }

  }
  #print(edgePot)

  # Now put log potentials into their corresponding places in a parameter vector:
  w <- numeric(crf$n.par)
  for(i in 1:nrow(nodePot)) {
    for(j in 1:ncol(nodePot)) {
      if(nodeMap[i,j,1] != 0) {
        w[nodeMap[i,j,1]] <- log(nodePot[i,j])
      }
    }
  }

  for(k in 1:length(edgePot)){
    for(i in 1:nrow(edgePot[[k]])){
      for(j in 1:ncol(edgePot[[k]])){
        if(edgeMap[[k]][i,j,1] != 0) {
          w[edgeMap[[k]][i,j,1]] <- log(edgePot[[k]][i,j])
        }
      }
    }
  }

  return(w)

}


#' Reformat parameter vector to CRF object forms node.pot and edge.pot BUT don't
#' exponentiate, i.e. leave them as logpots (energies).
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
par2logpots <- function(par, crf) {

  pots  <- make.pots(parms = par, crf = crf, rescaleQ = F, replaceQ = F, printQ = F)

  one.body.ens <- rep(list(NULL), nrow(pots[[1]]))
  for(i in 1:length(one.body.ens)){
    one.body.ens[[i]] <- log(pots[[1]])[i,,1]
  }

  two.body.ens <- rep(list(NULL), length(pots[[2]]))
  for(i in 1:length(pots[[2]])) {
    two.body.ens[[i]] <- log(pots[[2]][[i]])
  }

  return(list(one.body.ens, two.body.ens))
}
