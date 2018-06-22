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


#' Experiment: Modified Port of Schmidt UGM_MRF_NLL.m
#'
#' Experiment
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
loglik <- function(w, nInstances, suffStat, crf, inferFunc=infer.exact) {

  # Make potentials with the parameter vector passed in:
  mpots       <- make.pots(w, crf)
  crf$nodePot <- mpots$node.pot.shifted # Has to be done here because of arg structure to inferFunc
  crf$edgePot <- mpots$edge.pot.shifted # Has to be done here because of arg structure to inferFunc

  # Compute logZ:
  infer.info <- inferFunc(crf)
  logZZ      <- infer.info$logZ

  # Compute log likelihood:
  llk <- w %*% suffStat + nInstances*infer.info$logZ

  return(llk)
}


#' Utility function to compute gradient of log likelihood
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
grad.loglik <- function(crf, nInstances, suffStat, inference.func = infer.exact) {

  # First term of the gradient (its constant): N \times \hat{\text{E}}_{\boldsymbol \theta}[{\boldsymbol \phi}]
  #emp.num.features <- mrf.stat(crf, samples)
  emp.num.features <- suffStat

  # Second term of the gradient: N \times \text{E}_{\hat{\boldsymbol \theta}}[{\boldsymbol \phi}]
  #num.features.est <- nrow(samples) * feature.means(crf, inference.func)
  num.features.est <- nInstances * feature.means(crf, inference.func)

  grd <- emp.num.features - num.features.est

  return(grd)

}
