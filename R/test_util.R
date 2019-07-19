#' Joint distribution from true parameters for MRF
#'
#' Empirical joint distribution
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
fit_true <- function(true.mrf.modl) {

  potentials.info    <- make.gRbase.potentials(true.mrf.modl, node.names = 1:true.mrf.modl$n.nodes, state.nmes = c("1","2"))
  distribution.info  <- distribution.from.potentials(potentials.info$node.potentials, potentials.info$edge.potentials)
  joint.distribution <- as.data.frame(as.table(distribution.info$state.probs))
  #print(head(joint.distribution))

  # Re-order columns to increasing order
  freq.idx    <- ncol(joint.distribution)
  node.nums   <- colnames(joint.distribution)[-freq.idx]
  node.nums   <- unlist(strsplit(node.nums, split = "X"))
  node.nums   <- node.nums[-which(node.nums == "")]
  node.nums   <- as.numeric(node.nums)
  col.reorder <- order(node.nums)
  joint.distribution <- joint.distribution[,c(col.reorder, freq.idx)]

  return(joint.distribution)


}


#' Empirical joint distribution
#'
#' Empirical joint distribution
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
fit_empirical <- function(samples) {

  # Take a look at the empirical frequency-counts of the states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  config.counts      <- configs.and.counts[,freq.idx]

  # Estimate configuration state probabilities is by the EMPIRICAL relative frquencies:
  config.freqs                           <- config.counts/sum(config.counts)
  configs.and.counts                     <- cbind(configs.and.counts[-freq.idx],config.freqs)
  colnames(configs.and.counts)[freq.idx] <- "Emp.Freqs"

  # Order nodes in increasing order (in case they got out of wack):
  # NOTE: these should only be numbers for testing purposes!
  # BUT data.frame puts an X in front. Remove it to order
  node.nums   <- colnames(configs.and.counts)[-freq.idx]
  node.nums   <- unlist(strsplit(node.nums, split = "X"))
  node.nums   <- node.nums[-which(node.nums == "")]
  node.nums   <- as.numeric(node.nums)
  col.reorder <- order(node.nums)
  configs.and.counts <- configs.and.counts[,c(col.reorder, freq.idx)]

  return(configs.and.counts)

}


#' Joint distribution from MLE fit of parameters for full MRF loglikelihood
#'
#' Empirical joint distribution
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
fit_mle <- function(graph.eq, samples, inference.method = infer.exact, num.iter=10, mag.grad.tol=1e-3) {

  # Instantiate an empty model
  mle.mdl.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")

  # First compute the sufficient stats needed by the likelihood and its’ grad
  mle.mdl.fit$par.stat <- mrf.stat(mle.mdl.fit, samples)

  # Auxiliary, gradient convenience function. Follows train.mrf in CRF:
  gradient <- function(par, crf, ...) { crf$gradient }

  infr.meth <- inference.method     # inference method needed for Z and marginals calcs

  # May need to run the optimization a few times to get the gradient down:
  for(i in 1:num.iter) {
    opt.info  <- stats::optim(        # optimize parameters
      par          = mle.mdl.fit$par, # theta
      fn           = negloglik,       # objective function
      gr           = gradient,        # grad of obj func
      crf          = mle.mdl.fit,     # passed to fn/gr
      samples      = samples,         # passed to fn/gr
      infer.method = infr.meth,       # passed to fn/gr
      update.crfQ  = TRUE,            # passed to fn/gr
      method       = "L-BFGS-B",
      control      = list(trace = 1, REPORT=1))

    # Magnitude of the gradient:
    mag.grad <- sqrt(sum(mle.mdl.fit$gradient^2))

    print("==============================")
    print(paste("iter:", i))
    print(opt.info$convergence)
    print(opt.info$message)
    print(paste("Mag. Grad.:", mag.grad))
    print("==============================")

    if(mag.grad <= mag.grad.tol){
      break()
    }
  }
  #print(mle.mdl.fit$par)

  potentials.info    <- make.gRbase.potentials(mle.mdl.fit, node.names = colnames(samples), state.nmes = c("1","2"))
  distribution.info  <- distribution.from.potentials(potentials.info$node.potentials, potentials.info$edge.potentials)
  joint.distribution <- as.data.frame(as.table(distribution.info$state.probs))

  # Re-order columns to increasing order
  freq.idx    <- ncol(joint.distribution)
  node.nums   <- colnames(joint.distribution)[-freq.idx]
  node.nums   <- unlist(strsplit(node.nums, split = "X"))
  node.nums   <- node.nums[-which(node.nums == "")]
  node.nums   <- as.numeric(node.nums)
  col.reorder <- order(node.nums)
  joint.distribution <- joint.distribution[,c(col.reorder, freq.idx)]

  return(joint.distribution)

}
