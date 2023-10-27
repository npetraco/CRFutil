#' Compute the delta-alpha matrix {\boldsymbol \Delta}_{X_i=1}
#'
#' Needed for determination of parameters via logistic regression schemes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
delta.alpha <- function(crf, samples, state.nmes=c(1,2), printQ=FALSE, ff) {  # **** NEEDS TO BE C

  Da.mat        <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))
  cond.alp.mat  <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))
  cond.alpc.mat <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))

  count <- 1
  for(i in 1:crf$n.nodes) {
    for(n in 1:nrow(samples)) {

      X.cfg     <- samples[n,]
      Xc.cfg    <- X.cfg
      X.cfg[i]  <- state.nmes[1]  # Default state names are our usual choice.
      Xc.cfg[i] <- state.nmes[2]

      cond.alp  <- symbolic.conditional.energy(config = X.cfg,  condition.element.number = i, crf = crf, ff = ff, printQ = F, format = "conditional.phi")
      cond.alpc <- symbolic.conditional.energy(config = Xc.cfg, condition.element.number = i, crf = crf, ff = ff, printQ = F, format = "conditional.phi")
      cond.alp.mat[count,]  <- cond.alp
      cond.alpc.mat[count,] <- cond.alpc
      Da.mat[count,]      <- cond.alp - cond.alpc

      if(printQ == TRUE) {
        print(paste("Sample config#:", n, "Node:", i) )
        print("Configuration:")
        print(X.cfg)
        print("Complement Configuration:")
        print(Xc.cfg)
        print("Conditional alpha:")
        print(cond.alp)
        print("Complement Conditional alpha:")
        print(cond.alpc)
        print(paste("Delta-alpha", i, "Sample#:", n) )
        print(Da.mat[count,])
        print("=====================")
      }
      count <- count + 1
    }
  }


  da.info <- list(
    Da.mat,
    cond.alp.mat,
    cond.alpc.mat
  )
  names(da.info) <- c("Delta.alpha","conditional.alpha","comp.condtional.alpha")

  return(da.info)

}

#' Process Stan output for logistic fit
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
process.logistic.fit.stan <- function(stan.obj, param.vec.name, hpdi.unc=0.95) {

  par.mat  <- extract(stan.obj, param.vec.name)[[1]] # In package: rstan
  #num.pars <- ncol(par.mat)

  par.means <- colMeans(par.mat)
  par.meds  <- apply(par.mat, MARGIN = 2,FUN = median)
  par.intervals <- HPDinterval(as.mcmc(par.mat), prob = hpdi.unc) # In package: coda
  #print(par.intervals)

  par.info <- list(
    par.means,
    par.meds,
    par.intervals
  )
  names(par.info) <- c("par.post.means", "par.post.medians", "par.HPDIs")
  return(par.info)

}


#' Process JAGS output for logistic fit
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
process.logistic.fit.JAGS <- function(jags.obj) {

}
