#' Fit a MRF parameterized in the standard way
#'
#' Fit a MRF parameterized in the standard way
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


#' Fit a MRF via cross-validation glmnet LASSO logistic regression
#'
#' Fit a MRF via cross-validation glmnet LASSO logistic regression
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.glmnet.logistic.fit <- function(samples, crf.obj, lambda="lambda.1se", ff, infer.func, plotQ=F) {

  # Compute Delta-alpha matrix:
  print("Delta-alpha is too slow! FIX!!!!!!!!! See CRFutilRCppComponents")
  Delta.alpha.info.loc <- delta.alpha(crf = crf.obj, samples = as.matrix(samples), ff = ff)
  Delta.alpha.loc      <- Delta.alpha.info.loc$Delta.alpha
  print("Done finally........")

  # Prepare response vector.
  # ****NOTE: assumes states in the sample configurations are labeled 1 or 2!:
  y.loc <- as.numeric(samples)
  y.loc[which(y.loc==2)] <- 0

  # Logistic regression fit:
  glmnl.cv.loc <- glmnet::cv.glmnet(Delta.alpha.loc, y.loc, family = "binomial", type.measure = "class", intercept=F)

  # Process fit
  glmnl.fit.info <- glmnet2graph_info(glmnl.cv.loc, num.nodes = ncol(samples), lambda = lambda, plotQ = F)
  glmnl.edge.mat <- glmnl.fit.info$edge.parameters.mat[,c(3,4)]
  #print(glmnl.edge.mat)

  glmnl.adj      <- edges2adj(glmnl.edge.mat, n.nodes = ncol(samples), plotQ = F)
  #print(glmnl.adj)

  # Make another param est optional: none, mle, glmnet
  glmnl.crf      <- make.empty.field(adj.mat = glmnl.adj, parameterization.typ = "standard", plotQ = F)
  glmnl.out.pots <- make.pots(parms = glmnl.fit.info$parameter.vector, crf = glmnl.crf, replaceQ = T)

  #dump.crf(glmnl.crf)
  if(plotQ==T) {
    plot_crf(glmnl.crf)
  }

  return(glmnl.crf)

}
