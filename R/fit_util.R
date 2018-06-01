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
