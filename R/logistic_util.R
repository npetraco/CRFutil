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
delta.alpha <- function(crf, samples, printQ=FALSE) {

  Da.mat        <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))
  cond.alp.mat  <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))
  cond.alpc.mat <- array(NA,c(nrow(samples)*crf$n.nodes, crf$n.par))

  count <- 1
  for(i in 1:psl$n.nodes) {
    for(n in 1:nrow(samps)) {

      X.cfg     <- samples[n,]
      Xc.cfg    <- X.cfg
      X.cfg[i]  <- 1          # **** CAUTION: assumes states are labeled 1 or 2. GENERALIZE!!!!
      Xc.cfg[i] <- 2

      cond.alp  <- symbolic.conditional.energy(config = X.cfg,  condition.element.number = i, crf = crf, ff = f0, printQ = F, format = "conditional.phi")
      cond.alpc <- symbolic.conditional.energy(config = Xc.cfg, condition.element.number = i, crf = crf, ff = f0, printQ = F, format = "conditional.phi")
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
