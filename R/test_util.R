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


#' MLE fit of parameters ONLY for full MRF loglikelihood
#'
#' Use this for bigger models where it is not possible to compute the full joint distribution
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
fit_mle_params <- function(graph.eq, samples, parameterization.typ = "standard", opt.method="L-BFGS-B", inference.method = infer.exact, state.nmes = c("1","2"), num.iter=10, mag.grad.tol=1e-3, plotQ=F) {

  # Instantiate an empty model
  mle.mdl.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = parameterization.typ)

  # First compute the sufficient stats needed by the likelihood and its’ grad
  mle.mdl.fit$par.stat <- mrf.stat(mle.mdl.fit, samples)

  # Auxiliary, gradient convenience function. Follows train.mrf in CRF:
  gradient <- function(par, crf, ...) { crf$gradient }

  infr.meth <- inference.method     # inference method needed for Z and marginals calcs

  # May need to run the optimization a few times to get the gradient down:
  mag.grads <- NULL
  for(i in 1:num.iter) {
    opt.info  <- stats::optim(        # optimize parameters
      par          = mle.mdl.fit$par, # theta
      fn           = negloglik,       # objective function
      gr           = gradient,        # grad of obj func
      crf          = mle.mdl.fit,     # passed to fn/gr
      samples      = samples,         # passed to fn/gr
      infer.method = infr.meth,       # passed to fn/gr
      update.crfQ  = TRUE,            # passed to fn/gr
      method       = opt.method,
      control      = list(trace = 1, REPORT=1))

    # Magnitude of the gradient:
    mag.grad  <- sqrt(sum(mle.mdl.fit$gradient^2))
    mag.grads <- c(mag.grads, mag.grad)

    print("==============================")
    print(paste("iter:", i))
    print(opt.info$convergence)
    print(opt.info$message)
    print(paste("Mag. Grad.:", mag.grad))

    if(mag.grad <= mag.grad.tol){
      print("Gradient converged!")
      break()
    } else {
      print("Gradient not yet converged to specified tolerance")
    }
    print("==============================")

    if(plotQ==TRUE){
      plot(1:length(mag.grads), mag.grads, xlab="iteration", xlim=c(1,num.iter))
    }

  }
  #print(mle.mdl.fit$par)

  # Dress the potentials with gRbase decorations and include with crf object returned:
  out.pots        <- make.pots(parms = mle.mdl.fit$par,  crf = mle.mdl.fit,  rescaleQ = F, replaceQ = T)
  potentials.info <- make.gRbase.potentials(mle.mdl.fit, node.names = colnames(samples), state.nmes = state.nmes)
  mle.mdl.fit$node.potentials <- potentials.info$node.potentials
  mle.mdl.fit$edge.potentials <- potentials.info$edge.potentials
  mle.mdl.fit$node.energies   <- potentials.info$node.energies
  mle.mdl.fit$edge.energies   <- potentials.info$edge.energies

  #print(colnames(samples))

  return(mle.mdl.fit)

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


#' Joint distribution from bayes logistic regression fit of parameters
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
  logis.fit$par <- apply(extract(loc.bfit,"theta")[[1]], 2, median)
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
    loc.tmp                 <- node.potentials[[i]]
    loc.tmp                 <- loc.tmp/loc.tmp[2]               # Shift to standard parameterization
    loglin.fit$node.pot[i,] <- loc.tmp                      # Stick into mrf
    loglin.fit$par[count]   <- log(loglin.fit$node.pot[i,1]) # Stick into paramater vector
    count <- count + 1
  }

  # put edge potentials into mrf
  for(i in 1:length(edge.potentials)) {
    loc.tmp                  <- edge.potentials[[i]]
    dimnames(loc.tmp)        <- NULL                          # Strip dim names
    loc.tmp                  <- loc.tmp/loc.tmp[1,2]                  # Shift to standard parameterization
    loglin.fit$edge.pot[[i]] <- loc.tmp                           # Stick into mrf
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


#' Joint distribution from bayes poisson regression fit of parameters
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
fit_bayes_loglinear <- function(graph.eq, samples, iter=2000, thin=1, chains=4, control=NULL) {

  # Get frequency-counts of the configuration states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  loc.freqs          <- configs.and.counts[,freq.idx]
  loc.configs        <- configs.and.counts[,-freq.idx]

  # Construct model matrix. Try contr.sum version

  # Convert elements of all configs into factors. Required for model.matrix
  print("Building model matrix")
  loc.factor.mat <- apply(loc.configs,2, as.character)
  loc.factor.mat <- data.frame(loc.factor.mat)

  # Change contrasts:
  for(i in 1:ncol(loc.factor.mat)){
    contrasts(loc.factor.mat[,i]) <- contr.sum(2)
  }

  #model.matrix requires non numeric terms, so add Xs in formula
  loc.gfx <- adj2formula(ug(graph.eq, result = "matrix"), Xoption = T)
  loc.M   <- model.matrix(loc.gfx, data = loc.factor.mat)

  # Drop the intercept column. We handle it in Stan
  loc.M <- loc.M[,-1]
  #print(colnames(loc.M))
  #print(dim(loc.M))

  # Get frequency-counts of the configuration states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  loc.freqs          <- configs.and.counts[,freq.idx]

  loc.dat <- list(
    p = ncol(loc.M),
    N = nrow(loc.M),
    y = loc.freqs,
    Mmodl = loc.M
  )

  print("Compiling model")
  loc.model.c <- stanc(file = "inst/poisson_model.stan", model_name = 'model')
  loc.sm      <- stan_model(stanc_ret = loc.model.c, verbose = T)

  print("Sampling")
  loc.bfit <- sampling(loc.sm,
                       data    = loc.dat,
                       control = control,
                       iter    = iter,
                       thin    = thin,
                       chains  = chains)
  print("Done Sampling")

  # Instantiate an empty model
  bloglin.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")

  # Put coefs into mrf
  bloglin.fit$par <- apply(extract(loc.bfit,"theta")[[1]], 2, median)
  out.potsx       <- make.pots(parms = bloglin.fit$par, crf = bloglin.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(bloglin.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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


#' Joint distribution from bayes poisson regression fit of parameters, use phi model matrix
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
fit_bayes_loglinear2 <- function(graph.eq, samples, iter=2000, thin=1, chains=4, control=NULL, stan.model=NULL) {

  # Get frequency-counts of the configuration states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  loc.freqs          <- configs.and.counts[,freq.idx]
  loc.configs        <- configs.and.counts[,-freq.idx]

  # Instantiate an empty model
  bloglin.fit <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")
  loc.f0      <- function(yy){ as.numeric(c((yy==1),(yy==2)))}

  # Construct MRF model matrix.
  print("Building model matrix")
  loc.M  <- compute.model.matrix(
    configs   = loc.configs,
    edges.mat = bloglin.fit$edges,
    node.par  = bloglin.fit$node.par,
    edge.par  = bloglin.fit$edge.par,
    ff        = loc.f0)
  print("Done with model matrix. Sorry it's slow...")


  loc.dat <- list(
    p = ncol(loc.M),
    N = nrow(loc.M),
    y = loc.freqs,
    Mmodl = loc.M
  )

  if(is.null(stan.model)) {
    print("Compiling model")
    loc.model.c <- stanc(file = "inst/poisson_model.stan", model_name = 'model')
    loc.sm      <- stan_model(stanc_ret = loc.model.c, verbose = T)
  } else {
    loc.sm <- stan.model
  }

  print("Sampling")
  loc.bfit <- sampling(loc.sm,
                       data    = loc.dat,
                       control = control,
                       iter    = iter,
                       thin    = thin,
                       chains  = chains)
  print("Done Sampling")

  # Put coefs into mrf
  bloglin.fit$par <- apply(extract(loc.bfit,"theta")[[1]], 2, median)
  out.potsx       <- make.pots(parms = bloglin.fit$par, crf = bloglin.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(bloglin.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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


#' Joint distribution from bayes poisson regression fit of parameters, use phi model matrix
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
fit_bayes_zip <- function(graph.eq, samples, iter=2000, thin=1, chains=4, control=NULL, stan.model=NULL) {

  # Get frequency-counts of the configuration states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  loc.freqs          <- configs.and.counts[,freq.idx]
  loc.configs        <- configs.and.counts[,-freq.idx]

  # Instantiate an empty model
  bzipp.fit   <- make.empty.field(graph.eq = graph.eq, parameterization.typ = "standard")
  loc.f0      <- function(yy){ as.numeric(c((yy==1),(yy==2)))}

  # Construct MRF model matrix.
  print("Building model matrix")
  loc.M  <- compute.model.matrix(
    configs   = loc.configs,
    edges.mat = bzipp.fit$edges,
    node.par  = bzipp.fit$node.par,
    edge.par  = bzipp.fit$edge.par,
    ff        = loc.f0)
  print("Done with model matrix. Sorry it's slow...")


  # loc.dat <- list(
  #   K = ncol(loc.M),
  #   N = nrow(loc.M),
  #   y = loc.freqs,
  #   x = loc.M
  # )
  loc.dat <- list(
    M = ncol(loc.M),
    N = nrow(loc.M),
    y = loc.freqs,
    X = loc.M,
    s = rep(10,ncol(loc.M)),
    s_theta = rep(10,ncol(loc.M))
  )
  #print(loc.dat)

  if(is.null(stan.model)) {
    print("Compiling model")
    loc.model.c <- stanc(file = "inst/zero_inflated_poisson_take2.stan", model_name = 'model')
    loc.sm      <- stan_model(stanc_ret = loc.model.c, verbose = T)
  } else {
    loc.sm <- stan.model
  }

  print("Sampling")
  loc.bfit <- sampling(loc.sm,
                       data    = loc.dat,
                       control = control,
                       iter    = iter,
                       thin    = thin,
                       chains  = chains)
  print("Done Sampling")

  # Put coefs into mrf
  bzipp.fit$par   <- apply(extract(loc.bfit,"beta")[[1]], 2, median)
  # print("got here")

  out.potsx       <- make.pots(parms = bzipp.fit$par, crf = bzipp.fit, rescaleQ = T, replaceQ = T)

  potentials.info    <- make.gRbase.potentials(bzipp.fit, node.names = colnames(samples), state.nmes = c("1","2"))
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

  outinfo <- list(
    joint.distribution,
    loc.bfit
  )

  return(loc.bfit)

}


#' Generate a little lizard data model from Brunner course for testing
#'
#' XXXX
#'
#' The function will
#'
#' @param XX The XX
#'
#' @details The function will generate a little lizard data model for testing. Originally from: https://www.utstat.toronto.edu/~brunner/oldclass/312f12/lectures/312f12LoglinearWithR2.pdf
#'          Generating the data this way saves a little space in testing and allows the user
#'          to change the names of the nodes and states.
#'
#' @return A little data set
#'
#'
#' @export
generate_brunner_lizards <- function(typ=1, node1.name=NULL, node2.name=NULL, node3.name=NULL, node1.states=NULL, node2.states=NULL, node3.states=NULL) {

  lizards      <- numeric(8)
  dim(lizards) <- c(2,2,2) # Now a 2x2x2 table: Rows, cols, layers

  # 1 = Perch Height, 2 = Perch Diameter, 3 = Species
  lizards[,,1] <- rbind( c(15,18),
                         c(48,84) )
  lizards[,,2] <- rbind( c(21,1),
                         c(3, 2) )

  if(typ==1) {
    # Original labels Brunner used
    lizlabels          <- list()
    lizlabels$Height   <- c("gt_5.0","le_5.0")
    lizlabels$Diameter <- c("le_2.5","gt_2.5")
    lizlabels$Species  <- c("Sagrei","Angusticeps")
    dimnames(lizards)  <- lizlabels
  } else if(typ==2) {
    # More generic labels
    lizlabels     <- list()
    lizlabels$`1` <- c("1","2")
    lizlabels$`2` <- c("1","2")
    lizlabels$`3` <- c("1","2")
    dimnames(lizards)  <- lizlabels
  } else if(typ==3) {
    # User defined labels
    lizlabels     <- list(
      node1.states,
      node2.states,
      node3.states
    )
    names(lizlabels)  <- c(node1.name, node2.name, node3.name)
    dimnames(lizards) <- lizlabels
  } else if(typ==0) {
    # Do nothing for labels
    lizlabels     <- list()
  } else {
    stop("typ must == 0, 1, 2, or 3")
  }

  # Table forms
  lizards.config.freq.frame <- as.data.frame(as.table(lizards)) # Config freq table
  #print(lizards.config.freq.frame)

  # Configurations. Put in random order
  num.samp.configs <- sum(lizards)
  num.nodes.loc    <- (ncol(lizards.config.freq.frame)-1)
  lizards.configs  <- NULL

  for(i in 1:nrow(lizards.config.freq.frame)) {
    a.config    <- lizards.config.freq.frame[i, 1:num.nodes.loc]
    config.freq <- lizards.config.freq.frame[i, num.nodes.loc+1]

    for(j in 1:config.freq) {
      lizards.configs <- rbind(lizards.configs, a.config)
    }

  }
  rownames(lizards.configs) <- NULL
  #print(configs)

  # Empirical dist
  lizards.config.probs <- lizards/sum(lizards)
  lizards.config.probs <- as.data.frame(as.table(lizards.config.probs))

  lizards.info.list <- list(lizards,
                            lizards.config.freq.frame,
                            lizards.configs,
                            lizards.config.probs)

  names(lizards.info.list) <- c("lizards.contingency.table",
                                "lizards.configuration.frequencies",
                                "lizards.samples",
                                "lizards.configuration.probabilities")

  return(lizards.info.list)

}


#' Generate a little lizard data model from the gRbase package for testing
#'
#' XXXX
#'
#' The function will
#'
#' @param XX The XX
#'
#' @details The function will generate a little lizard data model for testing. Accessed from the gRbase package and
#' processed/reconfigured for various testing exercises.
#'
#' @return A little data set
#'
#'
#' @export
generate_gRbase_lizards <- function() {

  # gRbase lizard data should already be loaded
  #data(lizard)

  # Empirical dist
  lizards.config.probs <- lizard/sum(lizard)
  lizards.config.probs <- as.data.frame(as.table(lizards.config.probs))

  lizards.info.list <- list(lizard,    # Already in gRbase
                            lizardAGG, # Already in gRbase
                            lizardRAW, # Already in gRbase
                            lizards.config.probs)

  # Use the same names as for the Brunner data set above
  names(lizards.info.list) <- c("lizards.contingency.table",
                                "lizards.configuration.frequencies",
                                "lizards.samples",
                                "lizards.configuration.probabilities")

  return(lizards.info.list)

}


#' Generate the original triangle model used in Notes
#'
#' Generate the original triangle model used in Notes
#'
#' The function will generate the original triangle model used in Notes. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate the original triangle model used in Notes along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_triangle_model <- function(num.samples = NULL, Temp=1, plot.sampleQ=F, plot.graphQ=F) {

  # A:B + A:C + B:C
  tri.model.loc <- make.empty.field(graph.eq = ~1:2 + 1:3 + 2:3, parameterization.typ = "general", plot.graphQ)
  #dump.crf(tri.model.loc)

  # node potentials:
  tauA <- c( 1,    -1.3)
  tauB <- c(-0.85, -2.4)
  tauC <- c(3.82,   1.4)

  # edge potentials:
  omegaAB <- rbind(
    c( 3.5, -1.4),
    c(-1.4,  2.5)
  )
  omegaBC <- rbind(
    c( 2.6, 0.4),
    c( 0.4, 2.5)
  )
  omegaAC <- rbind(
    c(-0.6,  1.2),
    c( 1.2, -0.6)
  )

  tri.model.loc$node.pot <- rbind(
    exp(tauA/Temp),
    exp(tauB/Temp),
    exp(tauC/Temp)
  )
  rownames(tri.model.loc$node.pot) <- NULL # Remove the row names

  tri.model.loc$edge.pot[[1]] <- exp(omegaAB/Temp)
  tri.model.loc$edge.pot[[2]] <- exp(omegaAC/Temp)
  tri.model.loc$edge.pot[[3]] <- exp(omegaBC/Temp)

  tri.model.loc$par <- c(as.numeric(t(tri.model.loc$node.pot)), unlist(tri.model.loc$edge.pot))

  # Compute exact joint distribution from potentials:
  tri.model.ext.probs <- compute.full.distribution(tri.model.loc)$joint.distribution
  #print(tri.model.ext.probs)

  tri.model.samps             <- NULL
  tri.model.contingency       <- NULL
  tri.model.config.freq.frame <- NULL
  tri.model.probs             <- NULL
  if(!is.null(num.samples)) {
    tri.model.samps <- sample.exact(tri.model.loc, size = num.samples)
    #print(tri.model.samps)
    colnames(tri.model.samps) <- c("X1", "X2", "X3")

    if(plot.sampleQ == T) {
      mrf.sample.plot(tri.model.samps)
    }

    # Samples to contingency table
    tri.model.contingency <- xtabs(~., data=as.data.frame(tri.model.samps))
    #print(tri.model.contingency)

    # Sample config freq table
    tri.model.config.freq.frame <- as.data.frame(as.table(tri.model.contingency))

    # Sample configuration probabilities
    tri.model.emp.probs <- tri.model.contingency/sum(tri.model.contingency)
    tri.model.emp.probs <- as.data.frame(as.table(tri.model.emp.probs))
    tri.model.ext.probs <- align.distributions(tri.model.emp.probs, list(tri.model.ext.probs))
    #print(data.frame(tri.model.emp.probs, tri.model.ext.probs))


  }

  tri.model.info <-   list(
    tri.model.loc,
    tri.model.contingency,
    tri.model.config.freq.frame,
    tri.model.samps,
    tri.model.emp.probs,
    tri.model.ext.probs
  )

  names(tri.model.info) <- c("triangle.model",
                             "tri.contingency.table",
                             "tri.configuration.frequencies",
                             "tri.samples",
                             "tri.empirical.joint.distribution",
                             "tri.exact.joint.distribution")

  return(tri.model.info)

}


#' Generate the original Schmidt small chain model used in Notes
#'
#' Generate the original Schmidt small chain model used in Notes
#'
#' The function will generate the original Schmidt small chain model used in Notes. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate the original Schmidt small chain model used in Notes along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_schmidt_small_model <- function(num.samples = NULL, Temp = 1, plot.sampleQ=F, plot.graphQ=F) {

  # Cathy-Heather-Mark-Allison: 1-2-3-4
  model.loc <- make.empty.field(graph.eq = ~1:2 + 2:3 + 3:4, parameterization.typ = "general", plot.graphQ)
  #dump.crf(model.loc)

  Psi1 <- exp(log(c(0.25, 0.75)*4)/Temp)
  Psi2 <- exp(log(c(0.9,  0.1) *10)/Temp)
  Psi3 <- exp(log(c(0.25, 0.75)*4)/Temp)
  Psi4 <- exp(log(c(0.9,  0.1) *10)/Temp)

  Psi12 <- exp(log(6*rbind(c(2/6, 1/6),
                           c(1/6, 2/6)) )/Temp)
  Psi23 <- exp(log(6*rbind(c(2/6, 1/6),
                           c(1/6, 2/6)) )/Temp)
  Psi34 <- exp(log(6*rbind(c(2/6, 1/6),
                           c(1/6, 2/6)) )/Temp)

  model.loc$node.pot <- rbind(
    Psi1,
    Psi2,
    Psi3,
    Psi4
  )
  rownames(model.loc$node.pot) <- NULL # Remove the row names

  model.loc$edge.pot[[1]] <- Psi12
  model.loc$edge.pot[[2]] <- Psi23
  model.loc$edge.pot[[3]] <- Psi34

  model.loc$par <- c(as.numeric(t(model.loc$node.pot)), unlist(model.loc$edge.pot))
  #dump.crf(model.loc)

  # Compute exact joint distribution from potentials:
  model.ext.probs <- compute.full.distribution(model.loc)$joint.distribution
  #print(model.ext.probs)

  model.samps             <- NULL
  model.contingency       <- NULL
  model.config.freq.frame <- NULL
  model.probs             <- NULL
  if(!is.null(num.samples)) {
    model.samps <- sample.exact(model.loc, size = num.samples)
    colnames(model.samps) <- paste0("X", 1:model.loc$n.nodes)
    #print(head(model.samps))

    if(plot.sampleQ == T) {
      mrf.sample.plot(model.samps)
    }

    # Samples to contingency table
    model.contingency <- xtabs(~., data=as.data.frame(model.samps))
    #print(model.contingency)

    # Sample config freq table
    model.config.freq.frame <- as.data.frame(as.table(model.contingency))
    #print(model.config.freq.frame)

    # Sample configuration probabilities
    model.emp.probs <- model.contingency/sum(model.contingency)
    model.emp.probs <- as.data.frame(as.table(model.emp.probs))
    model.ext.probs <- align.distributions(model.emp.probs, list(model.ext.probs))
    #print(data.frame(model.emp.probs, model.ext.probs))

  }

  model.info <- list(
    model.loc,
    model.contingency,
    model.config.freq.frame,
    model.samps,
    model.emp.probs,
    model.ext.probs
  )

  names(model.info) <- c("schmidt.small.model",
                             "schmidt.small.contingency.table",
                             "schmidt.small.configuration.frequencies",
                             "schmidt.small.samples",
                             "schmidt.small.empirical.joint.distribution",
                             "schmidt.small.exact.joint.distribution")

  return(model.info)

}


#' Generate the original Koller-Friedman Misconception model used in their book and the Notes
#'
#' Generate the original Koller-Friedman Misconception model used in their book and the Notes
#'
#' The function will generate the original Koller-Friedman Misconception model used in Notes. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate the original Schmidt small chain model used in Notes along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_koller_misconception_model <- function(num.samples = NULL, Temp = 1, plot.sampleQ=F, plot.graphQ=F) {

  # Alice---Bob
  #    |     |
  # Charles-Debbie
  model.loc <- make.empty.field(graph.eq = ~1:2 + 2:3 + 3:4 + 4:1, parameterization.typ = "general", plotQ = plot.graphQ)
  #dump.crf(model.loc)

  # Node pots in original Misconception example are all 1 so no need to initialize them.
  # Edge pots in the Misconception example:
  PsiAB <- exp(log( rbind(
    c(30, 5),
    c(1, 10)) )/Temp)

  PsiBC <- exp(log( rbind(
    c(100, 1),
    c(1, 100)) )/Temp)

  PsiCD <- exp(log( rbind(
    c(1, 100),
    c(100, 1)) )/Temp)

  PsiDA <- exp(log( rbind(
    c(100, 1),
    c(1, 100)) )/Temp)

  model.loc$edge.pot[[1]] <- PsiAB
  model.loc$edge.pot[[2]] <- PsiBC
  model.loc$edge.pot[[3]] <- PsiCD
  model.loc$edge.pot[[4]] <- PsiDA

  model.loc$par <- c(as.numeric(t(model.loc$node.pot)), unlist(model.loc$edge.pot))
  #dump.crf(model.loc)

  # Compute exact joint distribution from potentials:
  model.ext.probs <- compute.full.distribution(model.loc)$joint.distribution
  #print(model.ext.probs)

  model.samps             <- NULL
  model.contingency       <- NULL
  model.config.freq.frame <- NULL
  model.probs             <- NULL
  if(!is.null(num.samples)) {
    model.samps <- sample.exact(model.loc, size = num.samples)
    colnames(model.samps) <- paste0("X", 1:model.loc$n.nodes)
    #print(head(model.samps))

    if(plot.sampleQ == T) {
      mrf.sample.plot(model.samps)
    }

    # Samples to contingency table
    model.contingency <- xtabs(~., data=as.data.frame(model.samps))
    #print(model.contingency)

    # Sample config freq table
    model.config.freq.frame <- as.data.frame(as.table(model.contingency))
    #print(model.config.freq.frame)

    # Sample configuration probabilities
    model.emp.probs <- model.contingency/sum(model.contingency)
    model.emp.probs <- as.data.frame(as.table(model.emp.probs))
    model.ext.probs <- align.distributions(model.emp.probs, list(model.ext.probs))
    #print(data.frame(model.emp.probs, model.ext.probs))

  }

  model.info <- list(
    model.loc,
    model.contingency,
    model.config.freq.frame,
    model.samps,
    model.emp.probs,
    model.ext.probs
  )

  names(model.info) <- c("koller.misconception.model",
                         "koller.misconception.contingency.table",
                         "koller.misconception.configuration.frequencies",
                         "koller.misconception.samples",
                         "koller.misconception.empirical.joint.distribution",
                         "koller.misconception.exact.joint.distribution")

  return(model.info)

}


#' Generate the star model used in their book and the Notes
#'
#' Generate the star model used in their book and the Notes
#'
#' The function will generate the star model used in Notes. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate the original Schmidt small chain model used in Notes along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_star_model <- function(num.samples = NULL, Temp=1, seed=NULL, plot.sampleQ=F, plot.graphQ=F) {


  # Graph formula for Star field:
  adj.loc   <- ug(~1:2+1:3+1:4+1:5+2:3+2:4+2:5+3:4+3:5, result="matrix")

  # Permute the rows and columns so that the node names are in ascending order. This makes it easier to check the edge potentials in the crf object
  n.nms   <- as.numeric(colnames(adj.loc))
  adj.loc <- adj.loc[order(n.nms), order(n.nms)]
  #print(adj.loc)

  model.loc <- make.empty.field(adj.mat = adj.loc, parameterization.typ = "standard", plotQ = plot.graphQ)

  if(!is.null(seed)) {
    set.seed(seed)
  }
  model.loc$par <- runif(model.loc$n.par,-1.5,1.5)/Temp
  out.pot.loc   <- make.pots(parms = model.loc$par,  crf = model.loc,  rescaleQ = F, replaceQ = T)
  #dump.crf(model.loc)

  # Compute exact joint distribution from potentials:
  model.ext.probs <- compute.full.distribution(model.loc)$joint.distribution
  #print(model.ext.probs)

  model.samps             <- NULL
  model.contingency       <- NULL
  model.config.freq.frame <- NULL
  model.probs             <- NULL
  if(!is.null(num.samples)) {
    model.samps <- sample.junction(model.loc, size = num.samples)
    colnames(model.samps) <- paste0("X", 1:model.loc$n.nodes)
    #print(head(model.samps))

    if(plot.sampleQ == T) {
      mrf.sample.plot(model.samps)
    }

    # Samples to contingency table
    model.contingency <- xtabs(~., data=as.data.frame(model.samps))
    #print(model.contingency)

    # Sample config freq table
    model.config.freq.frame <- as.data.frame(as.table(model.contingency))
    #print(model.config.freq.frame)

    # Sample configuration probabilities
    model.emp.probs <- model.contingency/sum(model.contingency)
    model.emp.probs <- as.data.frame(as.table(model.emp.probs))
    model.ext.probs <- align.distributions(model.emp.probs, list(model.ext.probs))
    #print(data.frame(model.emp.probs, model.ext.probs))

  }

  model.info <- list(
    model.loc,
    model.contingency,
    model.config.freq.frame,
    model.samps,
    model.emp.probs,
    model.ext.probs
  )

  names(model.info) <- c("star.model",
                         "star.contingency.table",
                         "star.configuration.frequencies",
                         "star.samples",
                         "star.empirical.joint.distribution",
                         "star.exact.joint.distribution")

  return(model.info)

}


#' Generate the tesseract model used in their book and the Notes
#'
#' Generate the tesseract model used in their book and the Notes
#'
#' The function will generate the tesseract model used in Notes. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate the tesseract model used in Notes along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_tesseract_model <- function(num.samples = NULL, Temp=1, seed=NULL, plot.sampleQ=F, plot.graphQ=F) {

  # Tesseract field model:
  adj.loc <- ug(~1:2  + 1:4   + 1:5  + 1:13 +
                2:3   + 2:6   + 2:14 +
                3:4   + 3:7   + 3:15 +
                4:8   + 4:16  +
                5:6   + 5:8   + 5:9 +
                6:7   + 6:10  +
                7:8   + 7:11  +
                8:12  +
                9:10  + 9:12  + 9:13 +
                10:11 + 10:14 +
                11:12 + 11:15 +
                12:16 +
                13:14 + 13:16 +
                14:15 +
                15:16, result="matrix")

  # Permute the rows and columns so that the node names are in ascending order. This makes it easier to check the edge potentials in the crf object
  n.nms   <- as.numeric(colnames(adj.loc))
  adj.loc <- adj.loc[order(n.nms), order(n.nms)]
  #print(adj.loc)

  model.loc <- make.empty.field(adj.mat = adj.loc, parameterization.typ = "standard", plotQ = plot.graphQ)

  if(!is.null(seed)) {
    set.seed(seed)
  }
  model.loc$par <- runif(model.loc$n.par,-1.5,1.5)/Temp
  out.pot.loc   <- make.pots(parms = model.loc$par,  crf = model.loc,  rescaleQ = F, replaceQ = T)
  #dump.crf(model.loc)

  # Compute exact joint distribution from potentials:
  model.ext.probs <- compute.full.distribution(model.loc)$joint.distribution
  #print(model.ext.probs)

  model.samps             <- NULL
  model.contingency       <- NULL
  model.config.freq.frame <- NULL
  model.probs             <- NULL
  if(!is.null(num.samples)) {
    model.samps <- sample.exact(model.loc, size = num.samples)
    colnames(model.samps) <- paste0("X", 1:model.loc$n.nodes)
    #print(head(model.samps))

    if(plot.sampleQ == T) {
      mrf.sample.plot(model.samps)
    }

    # Samples to contingency table
    model.contingency <- xtabs(~., data=as.data.frame(model.samps))
    #print(model.contingency)

    # Sample config freq table
    model.config.freq.frame <- as.data.frame(as.table(model.contingency))
    #print(model.config.freq.frame)

    # Sample configuration probabilities
    model.emp.probs <- model.contingency/sum(model.contingency)
    model.emp.probs <- as.data.frame(as.table(model.emp.probs))
    model.ext.probs <- align.distributions(model.emp.probs, list(model.ext.probs))
    #print(data.frame(model.emp.probs, model.ext.probs))

  }

  model.info <- list(
    model.loc,
    model.contingency,
    model.config.freq.frame,
    model.samps,
    model.emp.probs,
    model.ext.probs
  )

  names(model.info) <- c("tesseract.model",
                         "tesseract.contingency.table",
                         "tesseract.configuration.frequencies",
                         "tesseract.samples",
                         "tesseract.empirical.joint.distribution",
                         "tesseract.exact.joint.distribution")

  return(model.info)

}


#' Generate a random model using Erdos Renyi game
#'
#' Generate a random model using Erdos Renyi game
#'
#' The function will generate a random model using the Erdos Renyi game. Added a "Temp" parameter for user
#' to adjust the potentials in bulk. Also added sampling from the model.
#'
#' @param num.samples  Optional number of samples specification
#' @param Temp         Adjustment parameter to bulk adjust the potentials
#' @param plot.sampleQ Optional plot marginal node samples
#'
#' @details The function will generate a random model along with a random sample
#'
#' @return A little data set
#'
#'
#' @export
generate_random_model <- function(num.samples = NULL, Temp=1, num.nodes=10, p.or.numedges=0.6, type="gnp", seed=NULL, plot.sampleQ=F, plot.graphQ=T) {

  # A random graph:
  if(!is.null(seed)) {
    set.seed(seed)
  }
  gp.loc <- erdos.renyi.game(n = num.nodes, p.or.m = p.or.numedges, type = type)

  adj.loc <- as_adjacency_matrix(gp.loc, type = "both", sparse = F, names = T)
  # Name/rename the nodes to guarantee that they are in ascending order. This makes it easier to check the edge potentials in the crf object
  n.nms             <- 1:num.nodes
  colnames(adj.loc) <- n.nms
  rownames(adj.loc) <- n.nms
  #print(adj.loc)

  # Regenerate graph from adjacency matrix to get the node names right:
  # gp.loc  <- graph_from_adjacency_matrix(adj.loc, mode = "undirected")
  # if(plot.graphQ == T) {
  #   plot(gp.loc)
  # }

  model.loc <- make.empty.field(adj.mat = adj.loc, parameterization.typ = "standard", plotQ = plot.graphQ)

  if(is.null(seed)) {
    set.seed(seed)
  }
  model.loc$par <- runif(model.loc$n.par,-1.5,1.5)/Temp
  out.pot.loc   <- make.pots(parms = model.loc$par,  crf = model.loc,  rescaleQ = F, replaceQ = T)
  #dump.crf(model.loc)

  # Compute exact joint distribution from potentials:
  model.ext.probs <- compute.full.distribution(model.loc)$joint.distribution
  #print(model.ext.probs)

  model.samps             <- NULL
  model.contingency       <- NULL
  model.config.freq.frame <- NULL
  model.probs             <- NULL
  if(!is.null(num.samples)) {
    model.samps <- sample.exact(model.loc, size = num.samples)
    colnames(model.samps) <- paste0("X", 1:model.loc$n.nodes)
    #print(head(model.samps))

    if(plot.sampleQ == T) {
      mrf.sample.plot(model.samps)
    }

    # Samples to contingency table
    model.contingency <- xtabs(~., data=as.data.frame(model.samps))
    #print(model.contingency)

    # Sample config freq table
    model.config.freq.frame <- as.data.frame(as.table(model.contingency))
    #print(model.config.freq.frame)

    # Sample configuration probabilities
    model.emp.probs <- model.contingency/sum(model.contingency)
    model.emp.probs <- as.data.frame(as.table(model.emp.probs))
    model.ext.probs <- align.distributions(model.emp.probs, list(model.ext.probs))
    #print(data.frame(model.emp.probs, model.ext.probs))

  }

  model.info <- list(
    model.loc,
    model.contingency,
    model.config.freq.frame,
    model.samps,
    model.emp.probs,
    model.ext.probs
  )

  names(model.info) <- c("random.model",
                         "random.contingency.table",
                         "random.configuration.frequencies",
                         "random.samples",
                         "random.empirical.joint.distribution",
                         "random.exact.joint.distribution")

  return(model.info)

}
