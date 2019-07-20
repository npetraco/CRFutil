#' Joint distribution from true parameters for MRF
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
#' XXXXXX
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


#' Joint distribution from logistic regression fit of parameters
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
fit_logistic <- function(graph.eq, samples) {

  # Instantiate an empty model
  logis.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")

  # Make Delta-alpha matrix **** DELTA_ALPHA IS SLOW. IMPROVE
  print("Computing Delta-alpha")
  Delta.alpha.info <- delta.alpha(crf = logis.fit, samples = samples, printQ = F)
  Delta.alpha      <- Delta.alpha.info$Delta.alpha
  print("Done Delta-alpha. Sorry it's slow...")

  # Response vector
  y <- stack(data.frame(samples))[,1]
  y[which(y==2)] <- 0

  logis.glm.info <- glm(y ~ Delta.alpha -1, family=binomial(link="logit"))
  #print(summary(logis.glm.info))

  # Put coefs into mrf
  logis.fit$par <- as.numeric(coef(logis.glm.info))
  out.potsx     <- make.pots(parms = logis.fit$par, crf = logis.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(logis.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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


#' Joint distribution from logistic regression fit of parameters
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
fit_bayes_logistic <- function(graph.eq, samples, iter=2000, thin=1, chains=4, control=NULL) {

  # Instantiate an empty model
  logis.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")

  # Make Delta-alpha matrix **** DELTA_ALPHA IS SLOW. IMPROVE
  print("Computing Delta-alpha")
  Delta.alpha.info <- delta.alpha(crf = logis.fit, samples = samples, printQ = F)
  Delta.alpha      <- Delta.alpha.info$Delta.alpha
  print("Done Delta-alpha. Sorry it's slow...")

  # Response vector
  y <- stack(data.frame(samples))[,1]
  y[which(y==2)] <- 0

  #logis.glm.info <- glm(y ~ Delta.alpha -1, family=binomial(link="logit"))
  loc.model.c <- stanc(file = "inst/logistic_model.stan", model_name = 'model')
  print("Compiling model")
  loc.sm      <- stan_model(stanc_ret = loc.model.c, verbose = T)

  loc.dat <- list(
    N=nrow(Delta.alpha),
    K=ncol(Delta.alpha),
    Delta_alpha=Delta.alpha,
    y=y
  )

  print("Sampling")
  loc.bfit <- sampling(loc.sm,
                       data    = loc.dat,
                       control = control,
                       iter    = iter,
                       thin    = thin,
                       chains  = chains)
  print("Done Sampling")

  # Put coefs into mrf
  logis.fit$par <- apply(extract(bfit,"theta")[[1]], 2, median)
  out.potsx     <- make.pots(parms = logis.fit$par, crf = logis.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(logis.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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


#' Joint distribution from Poisson regression fit of parameters
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
fit_loglinear <- function(graph.eq, samples) {

  # Instantiate an empty model
  loglin.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")

  contingency.tab <- xtabs(~., data=samples)
  loglin.info     <- loglm(graph.eq, data=contingency.tab); # loglm from MASS which uses loglin in base
  loglin.coefs    <- coef(loglin.info)

  # Disentangle the loglin coefs so we can feed them in as potentials to distribution.from.potentials
  coef.clss <- sapply(1:length(loglin.coefs), function(xx){class(loglin.coefs[[xx]])})
  node.idxs <- which(coef.clss == "numeric")[-1] # first is the intercept
  edge.idxs <- which(coef.clss == "matrix")

  # Edges of the coefs. They are probably in a different order than in the crf$edges
  loglin.edges <- t(sapply(1:length(edge.idxs), function(xx){as.numeric(strsplit(x = names(loglin.coefs)[edge.idxs][xx], split = ".",fixed=T)[[1]])}))

  # Rearrange edges to be in the crf model order
  edg.rearr.idxs <- sapply(1:nrow(loglin.fit$edges), function(xx){row.match(x = loglin.fit$edges[xx,], loglin.edges)})
  # print(cbind(
  #   loglin.fit$edges,
  #   loglin.edges[edg.rearr.idxs,],
  #   names(loglin.coefs)[edge.idxs[edg.rearr.idxs] ],
  #   edge.idxs[edg.rearr.idxs],
  #   names(loglin.coefs)[edge.idxs]
  # ))

  # node potentials
  node.potentials <- lapply(loglin.coefs[node.idxs], FUN=exp)
  #print(node.potentials)

  # edge potentials
  edge.potentials <- lapply(loglin.coefs[edge.idxs[edg.rearr.idxs]], FUN=exp)

  # print(cbind(
  #   loglin.fit$edges,
  #   names(edge.potentials)
  # ))

  # put node potentials into mrf
  count <- 1
  for(i in 1:length(node.potentials)) {
    tmp                     <- node.potentials[[i]]
    tmp                     <- tmp/tmp[2]               # Shift to standard parameterization
    loglin.fit$node.pot[i,] <- tmp                      # Stick into mrf
    loglin.fit$par[count]   <- log(loglin.fit$node.pot[i,1]) # Stick into paramater vector
    count <- count + 1
  }

  # put edge potentials into mrf
  for(i in 1:length(edge.potentials)) {
    tmp                      <- edge.potentials[[i]]
    dimnames(tmp)            <- NULL                          # Strip dim names
    tmp                      <- tmp/tmp[1,2]                  # Shift to standard parameterization
    loglin.fit$edge.pot[[i]] <- tmp                           # Stick into mrf
    loglin.fit$par[count]    <- log(loglin.fit$edge.pot[[i]][1,1]) # Stick into paramater vector
    count <- count + 1
  }
  #print(loglin.fit$node.pot)
  #print(loglin.fit$edge.pot)
  #print(loglin.fit$par)

  out.potsx     <- make.pots(parms = loglin.fit$par, crf = loglin.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(loglin.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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
