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


#' Utility function to compute negative log pseudo-likelihood. Can also request the gradient.
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
neglogpseudolik.config <- function(param = NULL, config, crf, ff, cond.en.form="feature.function", gradQ=FALSE) {

  # Use input theta if supplied:
  if(is.null(param)){
    theta.pars <- crf$par
  } else {
    theta.pars <- param
  }

  num.nodes                       <- crf$n.nodes
  conditional.energies            <- array(NA, num.nodes)
  complement.conditional.energies <- array(NA, num.nodes)
  conditional.logZs               <- array(NA, num.nodes)
  if(gradQ==TRUE) {
    en.grad.mat             <- array(NA, c(length(theta.pars), num.nodes)) # Call alpha.mat??
    en.grad.mat.c           <- array(NA, c(length(theta.pars), num.nodes)) # Call alpha.mat.c??
    Z.grad.mat              <- array(NA, c(length(theta.pars), num.nodes))
    logZ.grad.mat           <- array(NA, c(length(theta.pars), num.nodes)) # Call Ealpha.mat??
    grad.neg.log.pseudo.lik <- numeric(length(theta.pars))
  }

  if(cond.en.form=="feature.function"){
    cond.en.func <- conditional.config.energy
  } else {
    if(cond.en.form=="feature"){
      cond.en.func <- conditional.config.energy2
    } else {
      stop("Conditional energy formula not properly specified! Use key word feature.function or feature for cond.en.form arguement.")
    }
  }
  for(i in 1:num.nodes) {

    # E(Xi | X/Xi)
    conditional.energies[i] <-
      cond.en.func(
        theta.par                = theta.pars,
        config                   = config,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)

    # E(Xi' | X/Xi')
    config.c     <- complement.at.idx(configuration = config, complement.index = i)
    complement.conditional.energies[i] <-
      cond.en.func(
        theta.par                = theta.pars,
        config                   = config.c,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)

    # Using above conditional energies, compute corresponding conditional logZ
    conditional.logZs[i] <- logsumexp2(c(conditional.energies[i], complement.conditional.energies[i]))

    if(gradQ==TRUE){
      # Conditional energy gradients, stored as columnns. These are \alpha_{[~i]}
      en.grad.mat[,i]   <- conditional.energy.gradient(config = config,   condition.element.number = i, crf = crf, ff = ff)$conditional.grad
      en.grad.mat.c[,i] <- conditional.energy.gradient(config = config.c, condition.element.number = i, crf = crf, ff = ff)$conditional.grad

      # Conditional partition function gradients, stored as columns.
      Z.grad.mat[,i] <- (exp(conditional.energies[i]) * en.grad.mat[,i] +
                         exp(complement.conditional.energies[i]) * en.grad.mat.c[,i])

      # Conditional log partition function gradients, stored as columns.
      # These are \text{E}[{\boldsymbol \alpha}_{[~i]}] = \nabla_{\boldsymbol \theta} \log(Z_{X_i|{\bf X}/X_i}) = \textbf{E}_{X_i} [{\boldsymbol \alpha}_{[\sim i]}]
      # with components \frac{\partial}{\partial \theta_k}\log\Big( Z_{X_i|{\bf X}\slash X_i} \Big) = \frac{1}{Z_{X_i|{\bf X}\slash X_i}} \frac{\partial}{\partial \theta_k} Z_{X_i|{\bf X}\slash X_i}
      logZ.grad.mat[,i] <- (1/exp(conditional.logZs[i])) * Z.grad.mat[,i] # COMBINE WITH ABOVE??

      # Gradient of neg log psedolikelihood for a configuration:
      # \nabla_{\boldsymbol \theta}{\cal L}_{\text{PL}}= \sum_{i=1}^p {\boldsymbol \alpha}_{[ \sim i]}({\bf X}) - \text{E}_{X_i} [{\boldsymbol \alpha}_{[\sim i]}]
      # do with colSums outside loop instead?
      #print(i)
      #print(grad.neg.log.pseudo.lik)
      # Note: gradient of the NEGATIVE log pseudo likelihood
      grad.neg.log.pseudo.lik <- grad.neg.log.pseudo.lik + (logZ.grad.mat[,i] - en.grad.mat[,i])
    }

  }

  # Note: NEGATIVE log pseudo likelihood
  neg.log.pseudo.lik <- sum(conditional.logZs - conditional.energies)

  if(gradQ==TRUE){
    neg.log.pseudo.lik.info <- list(
      neg.log.pseudo.lik,
      grad.neg.log.pseudo.lik
    )
  } else {
    neg.log.pseudo.lik.info <- list(
      neg.log.pseudo.lik,
      NULL
    )
  }

  #----
  # intermediate.info <- list(
  #   conditional.logZs,
  #   en.grad.mat,   #  alpha.mat
  #   en.grad.mat.c, #  alpha.mat.c
  #   Z.grad.mat,
  #   logZ.grad.mat #  E[alpha.mat]
  # )
  # names(intermediate.info) <- c("logZs","alpha","alpha.c","dZ","Ealpha")
  #----

  names(neg.log.pseudo.lik.info) <- c("neglogpseudolik","grad.neglogpseudolik")

  # return(list(neg.log.pseudo.lik.info, intermediate.info))
  return(neg.log.pseudo.lik.info)

}


#' Utility function to SEPARATELY compute gradient negative log pseudo-likelihood.
#'
#' Assumes features are 0,1 valued and parameters are numbered.
#'
#' @details neglogpseudolik.config() computes the gradient as as well. This is just a pull out for testing and pedegogy.
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
grad.neglogpseudolik.config <- function(param = NULL, config, crf, ff, cond.en.form="feature.function") {

  # Use input theta if supplied:
  if(is.null(param)){
    theta.pars <- crf$par
  } else {
    theta.pars <- param
  }

  num.nodes                       <- crf$n.nodes
  conditional.energies            <- array(NA, num.nodes)
  complement.conditional.energies <- array(NA, num.nodes)
  conditional.logZs               <- array(NA, num.nodes)
  en.grad.mat             <- array(NA, c(length(theta.pars), num.nodes)) # Call alpha.mat??
  en.grad.mat.c           <- array(NA, c(length(theta.pars), num.nodes)) # Call alpha.mat.c??
  Z.grad.mat              <- array(NA, c(length(theta.pars), num.nodes))
  logZ.grad.mat           <- array(NA, c(length(theta.pars), num.nodes)) # Call Ealpha.mat??
  grad.neg.log.pseudo.lik <- numeric(length(theta.pars))

  for(i in 1:num.nodes) {

    if(cond.en.form=="feature.function"){
      cond.en.func <- conditional.config.energy
    } else {
      if(cond.en.form=="feature"){
        cond.en.func <- conditional.config.energy2
      } else {
        stop("Conditional energy formula not properly specified! Use key word feature.function or feature for cond.en.form arguement.")
      }
    }

    # E(Xi | X/Xi)
    conditional.energies[i] <-
      cond.en.func(
        theta.par                = theta.pars,
        config                   = config,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)

    # E(Xi' | X/Xi')
    config.c     <- complement.at.idx(configuration = config, complement.index = i)
    complement.conditional.energies[i] <-
      cond.en.func(
        theta.par                = theta.pars,
        config                   = config.c,
        condition.element.number = i,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)

    # Using above conditional energies, compute corresponding conditional logZ
    conditional.logZs[i] <- logsumexp2(c(conditional.energies[i], complement.conditional.energies[i]))

    # Conditional energy gradients, stored as columnns. These are \alpha_{[~i]}
    en.grad.mat[,i]   <- conditional.energy.gradient(config = config,   condition.element.number = i, crf = crf, ff = ff)$conditional.grad
    en.grad.mat.c[,i] <- conditional.energy.gradient(config = config.c, condition.element.number = i, crf = crf, ff = ff)$conditional.grad

    # Conditional partition function gradients, stored as columns.
    Z.grad.mat[,i] <- (exp(conditional.energies[i]) * en.grad.mat[,i] +
                         exp(complement.conditional.energies[i]) * en.grad.mat.c[,i])

    # Conditional log partition function gradients, stored as columns.
    # These are \text{E}[{\boldsymbol \alpha}_{[~i]}] = \nabla_{\boldsymbol \theta} \log(Z_{X_i|{\bf X}/X_i}) = \textbf{E}_{X_i} [{\boldsymbol \alpha}_{[\sim i]}]
    # with components \frac{\partial}{\partial \theta_k}\log\Big( Z_{X_i|{\bf X}\slash X_i} \Big) = \frac{1}{Z_{X_i|{\bf X}\slash X_i}} \frac{\partial}{\partial \theta_k} Z_{X_i|{\bf X}\slash X_i}
    logZ.grad.mat[,i] <- (1/exp(conditional.logZs[i])) * Z.grad.mat[,i] # COMBINE WITH ABOVE??

    # Gradient of neg log psedolikelihood for a configuration:
    # \nabla_{\boldsymbol \theta}{\cal L}_{\text{PL}}= \sum_{i=1}^p {\boldsymbol \alpha}_{[ \sim i]}({\bf X}) - \text{E}_{X_i} [{\boldsymbol \alpha}_{[\sim i]}]
    # do with colSums outside loop instead?
    #print(i)
    #print(grad.neg.log.pseudo.lik)
    # Note: gradient of the NEGATIVE log pseudo likelihood
    grad.neg.log.pseudo.lik <- grad.neg.log.pseudo.lik + (logZ.grad.mat[,i] - en.grad.mat[,i])
  }

  #----
  grad.info <- list(
    grad.neg.log.pseudo.lik, #
    conditional.logZs,
    en.grad.mat,             #  alpha.mat
    en.grad.mat.c,           #  alpha.mat.c
    Z.grad.mat,
    logZ.grad.mat            #  E[alpha.mat]
  )
  names(grad.info) <- c("grad.neglogpseudolik","logZs","alpha","alpha.c","dZ","Ealpha")
  #----

  return(grad.info)

}


#' Sample negative log pseudo-likelihood
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
neglogpseudolik <- function(param, crf, samples, conditional.energy.function.type="feature.function", ff, gradQ = FALSE, update.crfQ = TRUE) {

  # Use input theta if supplied:
  if(is.null(param)){
    theta.pars <- crf$par
  } else {
    theta.pars <- param
  }

  # **** Good candidate for a parallelized loop????
  nliks  <- array(NA,nrow(samples))
  gnliks <- array(NA,c(nrow(samples),length(theta.pars)))
  nlik <- 0
  gnlik <- numeric(length(theta.pars))
  for(i in 1:nrow(samples)) {

    #print(paste("Sample:",i))
    nlik.info  <- neglogpseudolik.config(param = theta.pars, config = samples[i,], crf = crf, ff = ff, cond.en.form = conditional.energy.function.type, gradQ = gradQ)
    nliks[i]   <- nlik.info$neglogpseudolik
    nlik       <- nlik + nliks[i]
    if(gradQ==TRUE){
      gnliks[i,] <- nlik.info$grad.neglogpseudolik
      gnlik <- gnlik + gnliks[i,]
    }
    #print(nliks[i])
    #print(gnliks[i,])

  }

  if(gradQ==TRUE){
    samp.neg.log.pseudo.lik.info <- list(
      nlik,
      gnlik,
      nliks,
      gnliks
    )
  } else {
    samp.neg.log.pseudo.lik.info <- list(
      nlik,
      NULL,
      nliks,
      NULL
    )
  }

  if(update.crfQ == TRUE){
    crf$par      <- theta.pars
    crf$nll      <- nlik
    crf$gradient <- gnlik
  }

  names(samp.neg.log.pseudo.lik.info) <- c("samp.neglogpseudolik","samp.grad.neglogpseudolik","nliks","gnliks")

  return(samp.neg.log.pseudo.lik.info)

}


#' Sample negative log pseudo-likelihood
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
neglogpseudolik2 <- function(param, crf, samples, conditional.energy.function.type="feature.function", ff, gradQ = FALSE, update.crfQ = TRUE) {

  # Use input theta if supplied:
  if(is.null(param)){
    theta.pars <- crf$par
  } else {
    theta.pars <- param
  }

  # **** Good candidate for a parallelized loop????
  nliks  <- array(NA,nrow(samples))
  gnliks <- array(NA,c(nrow(samples),length(theta.pars)))
  nlik <- 0
  gnlik <- numeric(length(theta.pars))
  for(i in 1:nrow(samples)) {

    #print(paste("Sample:",i))
    nlik.info  <- neglogpseudolik.config(param = theta.pars, config = samples[i,], crf = crf, ff = ff, cond.en.form = conditional.energy.function.type, gradQ = gradQ)
    nliks[i]   <- nlik.info$neglogpseudolik
    nlik       <- nlik + nliks[i]
    if(gradQ==TRUE){
      gnliks[i,] <- nlik.info$grad.neglogpseudolik
      gnlik <- gnlik + gnliks[i,]
    }
    #print(nliks[i])
    #print(gnliks[i,])

  }

  if(gradQ==TRUE){
    samp.neg.log.pseudo.lik.info <- list(
      nlik,
      gnlik,
      nliks,
      gnliks
    )
  } else {
    samp.neg.log.pseudo.lik.info <- list(
      nlik,
      NULL,
      nliks,
      NULL
    )
  }

  if(update.crfQ == TRUE){
    crf$par      <- theta.pars
    crf$nll      <- nlik
    crf$gradient <- gnlik
  }

  #names(samp.neg.log.pseudo.lik.info) <- c("samp.neglogpseudolik","samp.grad.neglogpseudolik","nliks","gnliks")

  return(nlik)

}


#' Gradient of sample negative log pseudo-likelihood
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
grad.neglogpseudolik.FIXME <- function(par, crf, samples, ff) {

  grd <- rowSums(sapply(1:nrow(samples), function(xx){
    grad.neglogpseudolik.config(
      config = samples[xx,],
      phi.config = NULL,
      node.conditional.energies = NULL,
      node.complement.conditional.energies = NULL,
      par = par,
      crf = crf,
      ff = ff)$gradient.log.pseudolikelihood
    }))

  return(grd)

}
