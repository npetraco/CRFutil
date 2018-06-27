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
negloglik <- function(par, crf, samples, infer.method = infer.exact, update.crfQ = TRUE) {
        # args were: (w, nInstances, suffStat, crf, inferFunc=infer.exact)
  # args to mrf.nll: (par, crf, instances, infer.method = infer.chain, ...)


  # Make potentials with the parameter vector passed in.
  # Store in crf object for use with infer.method
  # They are needed to get the Z corresponding to w passed in!
  mpots <- make.pots(par, crf, rescaleQ = F, replaceQ = T, printQ = F)

  # Compute logZ:
  infer.info <- infer.method(crf)
  logZZ      <- infer.info$logZ

  # Compute neg log likelihood:
  if(is.null(crf$par.stat)) {
    stop("Compute sufficient statistics and store in crf object!")
  }
  nInstances <- nrow(samples)
  suffStat   <- crf$par.stat
  nllk       <- -par %*% suffStat + nInstances * infer.info$logZ

  # Compute gradient of neg log likelihood, just like in mrf.XXXX.nll():
  grad <- grad.negloglik(crf, nInstances, suffStat, infer.method)

  if(update.crfQ == TRUE){
    crf$par      <- par
    crf$nll      <- nllk
    crf$gradient <- grad
  }

  return(nllk)
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
grad.negloglik <- function(crf, nInstances, suffStat, inference.func = infer.exact) {

  # Second term of the gradient (its constant): N \times \hat{\text{E}}_{\boldsymbol \theta}[{\boldsymbol \phi}]
  #emp.num.features <- mrf.stat(crf, samples)
  emp.num.features <- suffStat

  # First term of the gradient: N \times \text{E}_{\hat{\boldsymbol \theta}}[{\boldsymbol \phi}]
  #num.features.est <- nrow(samples) * feature.means(crf, inference.func)
  num.features.est <- nInstances * feature.means(crf, inference.func)

  grd <- num.features.est - emp.num.features

  return(grd)

}
