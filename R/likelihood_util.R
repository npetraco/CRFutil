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


#' Utility function to compute negative log pseudo-likelihood
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
neglogpseudolik.config <- function(par = NULL, config, crf, ff) {

  num.nodes                      <- crf$n.nodes
  condtional.energies            <- array(NA, num.nodes)
  complement.condtional.energies <- array(NA, num.nodes)
  conditional.logZs              <- array(NA, num.nodes)

  if(is.null(par)){
    theta.pars <- crf$par
  } else {
    theta.pars <- par
  }

  phi.config <- phi.features(
    config    = config,
    edges.mat = crf$edges,
    node.par  = crf$node.par,
    edge.par  = crf$edge.par,
    ff        = ff)

  for(i in 1:num.nodes) {

    condtional.energies[i] <- conditional.config.energy2a(
      par                      = theta.pars,
      config                   = config,
      phi.config               = phi.config,
      condition.element.number = i,
      crf                      = crf,
      ff                       = ff)

    config.c     <- complement.at.idx(configuration = config, complement.index = i)
    phi.config.c <- phi.features(
      config    = config.c,
      edges.mat = crf$edges,
      node.par  = crf$node.par,
      edge.par  = crf$edge.par,
      ff        = ff)

    complement.condtional.energies[i] <- conditional.config.energy2a(
      par                      = theta.pars,
      config                   = config.c,
      phi.config               = phi.config.c,
      condition.element.number = i,
      crf                      = crf,
      ff                       = ff)

    conditional.logZs[i] <- logsumexp2(c(condtional.energies[i], complement.condtional.energies[i]))

  }

  #neg.log.pseudo.lik <- -sum(condtional.energies - log(conditional.Zs))
  neg.log.pseudo.lik <- sum(conditional.logZs - condtional.energies)

  return(neg.log.pseudo.lik)

}


#' Utility function to compute gradient of negative log pseudo-likelihood
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
conditional.config.energy.grad.test <- function(config, crf, ff) {

  if(is.null(crf$nodes2pars)){
    stop("Compute node-parameter associations and store in crf object!")
  }

  config.phi.vec <- phi.features(
    config    = config,
    edges.mat = crf$edges,
    node.par  = crf$node.par,
    edge.par  = crf$edge.par,
    ff        = ff
  )

  # Gradient "matrix" is #parameters by #nodes. I.E., each column is a gradient of a condtional energy:
  gradient.mat <- array(NA, c(crf$n.par, crf$n.nodes))
  for(i in 1:crf$n.nodes) {

    node.pars             <- crf$nodes2pars[[i]]       # Definitley derivs NOT with respect to these params are 0
    dEconfig.i            <- numeric(crf$n.par)        # Initalize a conditional energy gradient vector for node i to 0s
    dEconfig.i[node.pars] <- config.phi.vec[node.pars] # Any phi_i = 0 in here are also 0 derivs

    gradient.mat[,i] <- dEconfig.i
  }
  colnames(gradient.mat) <- 1:crf$n.nodes
  rownames(gradient.mat) <- 1:crf$n.par

  return(gradient.mat)

}
