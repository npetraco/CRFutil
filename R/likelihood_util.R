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
  conditional.energies            <- array(NA, num.nodes)
  complement.conditional.energies <- array(NA, num.nodes)
  conditional.logZs              <- array(NA, num.nodes)
  #----
  # conditional.Zs                 <- array(NA, num.nodes)
  # conditional.Prs                <- array(NA, num.nodes)
  # complement.conditional.Prs     <- array(NA, num.nodes)
  #----

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

    conditional.energies[i] <- conditional.config.energy2a(
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

    complement.conditional.energies[i] <- conditional.config.energy2a(
      par                      = theta.pars,
      config                   = config.c,
      phi.config               = phi.config.c,
      condition.element.number = i,
      crf                      = crf,
      ff                       = ff)

    #----
    #conditional.Zs[i]    <- exp(conditional.energies[i]) + exp(complement.conditional.energies[i])
    #complement.conditional.Prs[i] <- 1 - conditional.Prs[i]
    #----

    conditional.logZs[i] <- logsumexp2(c(conditional.energies[i], complement.conditional.energies[i]))

  }

  #----
  # print("CE")
  # print(conditional.energies)
  # print("CCE")
  # print(complement.conditional.energies)
  #
  # print("Zc")
  # print(conditional.Zs)
  #
  # print("cPrs")
  # print(conditional.Prs)
  # print("ccPrs")
  # print(complement.conditional.Prs)
  #----

  #neg.log.pseudo.lik <- -sum(conditional.energies - log(conditional.Zs))
  neg.log.pseudo.lik <- sum(conditional.logZs - conditional.energies)

  #----
  #take2 <- -log(prod(conditional.Prs))
  # print(take2)
  #
  #take3 <- -sum(log(conditional.Prs))
  # print(take3)
  #take4 <- -sum(conditional.energies - log(conditional.Zs))
  # print(take4)
  #

  # intermediate.info <- list(
  #   conditional.energies,
  #   complement.conditional.energies,
  #   conditional.Zs,
  #   conditional.logZs,
  #   conditional.Prs,
  #   complement.conditional.Prs,
  #   neg.log.pseudo.lik,
  #   take2,take3,take4)
  #----

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
grad.neglogpseudolik.config <- function(config, phi.config=NULL, node.conditional.energies=NULL, node.complement.conditional.energies=NULL, par=NULL, crf, ff) {

  if(is.null(crf$nodes2pars)) {
    stop("Compute node to parameter associations and store in crf object!")
  }

  # Best to pass these in so that they are not needlessly re-computed if this function is in a loop.
  # However, if they're not passed in, we need to compute them.
  if(is.null(phi.config)) {
    phi.config <- phi.features(
      config    = config,
      edges.mat = crf$edges,
      node.par  = crf$node.par,
      edge.par  = crf$edge.par,
      ff        = f0
    )
  }
  if(is.null(node.conditional.energies)){
    conditional.energies <- numeric(crf$n.nodes)
  } else {
    conditional.energies <- node.conditional.energies
  }
  if(is.null(node.complement.conditional.energies)){
    complement.conditional.energies <- numeric(crf$n.nodes)
  } else {
    complement.conditional.energies <- node.complement.conditional.energies
  }

  # If no parameter vector is supplied, use the one in the crf object
  if(is.null(par)) {
    theta.par <- crf$par
  } else {
    theta.par <- par
  }

  # Storage containers for requisite intermediates:
  # Conditional Energy Gradient "matrices": #parameters by #nodes. I.E., each column is a gradient of a condtional energy:
  grad.cond.E.mat      <- array(NA, c(crf$n.par, crf$n.nodes))
  grad.comp.cond.E.mat <- array(NA, c(crf$n.par, crf$n.nodes))

  # Condtional and complement condtional energies: vector length = #nodes
  cond.E.vec      <- array(NA, c(crf$n.nodes))
  comp.cond.E.vec <- array(NA, c(crf$n.nodes))

  # Condtional partition functions: vector length = #nodes
  cond.Z.vec  <- array(NA, c(crf$n.nodes))

  # *********ELIMINATE???????
  # Gradients of condtional partition functions: (\nabla_{\boldsymbol \theta}\textit{\textbf{Z}})_{ki}
  # Conditional partition function gradient "matrices": #parameters by #nodes.
  # I.E., each column is a gradient of a conditional partition function:
  grad.cond.Z.mat  <- array(NA, c(crf$n.par, crf$n.nodes))

  # Gradients of log condtional partition functions: \Big(\nabla_{\boldsymbol \theta} \log\big( \textbf{\textit{Z}}_{X_i|{\bf X}\slash X_i} \big)\Big)_{ki} = \Big( \textbf{E}_{{\bf X}} [{\boldsymbol \phi}_{[\sim i]}] \Big)_{ki}
  E.phi.config.mat <- array(NA, c(crf$n.par, crf$n.nodes))

  for(i in 1:crf$n.nodes) {

    # ** IMPORTANT PART ** Definitley derivs NOT with respect to these params are 0:
    node.pars <- crf$nodes2pars[[i]]

    # Initalize conditional energy gradient vectors for node i to 0s:
    dEX.i  <- numeric(crf$n.par)
    dEXc.i <- numeric(crf$n.par)

    # Get complement phi_i's needed for conditional Z's and log Z's:
    config.c     <- complement.at.idx(config,i)
    phi.config.c <- phi.features(
      config    = config.c,
      edges.mat = crf$edges,
      node.par  = crf$node.par,
      edge.par  = crf$edge.par,
      ff        = ff
    )

    # ** IMPORTANT PART ** Any phi_i = 0 in here are also 0 derivs:
    dEX.i[node.pars]  <- phi.config[node.pars]
    dEXc.i[node.pars] <- phi.config.c[node.pars]

    # Store condional energy gradients column-wise:
    grad.cond.E.mat[,i]      <- dEX.i
    grad.comp.cond.E.mat[,i] <- dEXc.i

    # Compute conditional energies if they were not passed in:
    if(is.null(node.conditional.energies)){
      conditional.energies[i] <- conditional.config.energy2a(
        par                      = theta.par,
        config                   = config,
        phi.config               = phi.config,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff)
    }
    if(is.null(node.complement.conditional.energies)) {
      complement.conditional.energies[i] <- conditional.config.energy2a(
        par                      = theta.par,
        config                   = config.c,
        phi.config               = phi.config.c,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff)
    }

    # Store conditional Z gradients column-wise as well:
    # NOTE: Assumes only two states per node
    grad.cond.Z.mat[,i] <- exp(conditional.energies[i]) * grad.cond.E.mat[,i] + exp(complement.conditional.energies[i]) * grad.comp.cond.E.mat[,i]

    # Condtional Zs:
    # NOTE: Assumes only two states per node
    cond.Z.vec[i] <- exp(conditional.energies[i]) + exp(complement.conditional.energies[i])

    # Average features with respect to states of X_i: \textbf{E}_{{\bf X}} [{\boldsymbol \phi}_{[\sim i]}] = \nabla_{\boldsymbol \theta} \log\big( \textbf{\textit{Z}}_{X_i|{\bf X}\slash X_i} \big)
    E.phi.config.mat[,i] <- grad.cond.Z.mat[,i]/cond.Z.vec[i]

  }

  # Gradient of pseudolikelihood for a config: \nabla_{\boldsymbol \theta} {\cal L}_{\text{PL}}
  grad.psl.config <- rowSums(grad.cond.E.mat - E.phi.config.mat)

  colnames(grad.cond.E.mat)      <- 1:crf$n.nodes
  rownames(grad.cond.E.mat)      <- 1:crf$n.par
  colnames(grad.comp.cond.E.mat) <- 1:crf$n.nodes
  rownames(grad.comp.cond.E.mat) <- 1:crf$n.par
  colnames(grad.cond.Z.mat)      <- 1:crf$n.nodes
  rownames(grad.cond.Z.mat)      <- 1:crf$n.par
  colnames(E.phi.config.mat)     <- 1:crf$n.nodes
  rownames(E.phi.config.mat)     <- 1:crf$n.par

  grad.info <- list(
    conditional.energies,
    complement.conditional.energies,
    grad.cond.E.mat,
    grad.comp.cond.E.mat,
    grad.cond.Z.mat,
    E.phi.config.mat,
    grad.psl.config
  )

  names(grad.info) <- c(
    "conditional.energies",
    "complement.conditional.energies",
    "gradients.conditional.energies",
    "gradients.complement.conditional.energies",
    "gradients.conditional.partition.functions",
    "gradients.log.conditional.partition.functions",
    "gradient.log.pseudolikelihood"
  )

  return(grad.info)

  #return(grad.psl.config)

}


#' TEST Utility function to compute gradient of negative log pseudo-likelihood
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
