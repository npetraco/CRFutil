#' Fit a MRF paramiterized in the standard way
#'
#' Fit a MRF paramiterized in the standard way
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.standard.fit <- function(samples, net.graph.formula, num.states, mrf.nll.func, infer.func) {

  adj <- ug(net.graph.formula, result="matrix")

  mrf.fit <- make.crf(adj, num.states)
  mrf.fit <- make.features(mrf.fit)

  # Standard parameterization:
  num.params             <- mrf.fit$n.nodes + mrf.fit$n.edges
  mrf.fit                <- make.par(mrf.fit, num.params)
  mrf.fit$node.par[,1,1] <- 1:mrf.fit$n.nodes

  for(i in 1:mrf.fit$n.edges){
    mrf.fit$edge.par[[i]][1,1,1] <- mrf.fit$n.nodes + i
    mrf.fit$edge.par[[i]][2,2,1] <- mrf.fit$n.nodes + i
  }

  mrf.fit <- train.mrf(mrf.fit, nll = mrf.nll.func, samples, infer.method = infer.func)
  infer.info.fit <- infer.func(mrf.fit)
  logZZ <- infer.info.fit$logZ

  fit.info <- list(mrf.fit, infer.info.fit)

  names(fit.info) <- c("fit.model","inference.info")

  return(fit.info)

}


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
  print("Potentials shifted to parameter vector via node.par and edge.par")

}


#' Extract parameter vector theta from node and potentials written by trin.mrf
#'
#' For testing purposes. train.mrf par and potentials don't always match exactly because of shifting of the pots
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
extract.theta.from.pots <- function(crf) {

}
