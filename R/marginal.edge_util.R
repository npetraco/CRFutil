#' marginal.edge.mrf
#'
#' XXXXXXX
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.mrf <- function(edge.samples){

  if(ncol(edge.samples) != 2) {
    stop("Input one edge (two nodes) only!")
  }

  # Construct contingency table for edge data
  edge.samps.loc           <- edge.samples # Make a copy
  if(is.null(colnames(edge.samples))){
    colnames(edge.samps.loc) <- c("X1","X2")
  } else {
    colnames(edge.samps.loc) <- colnames(edge.samples)
  }

  # Fit a two node MRF to examine the edge potentials and resulting node/edge beliefs.
  edge.grphf         <- ~EX1:EX2
  edge.adj           <- ug(edge.grphf, result="matrix")
  nnmes              <- colnames(edge.samps.loc)
  colnames(edge.adj) <- nnmes
  rownames(edge.adj) <- nnmes
  edge.n.nodes       <- 2
  edge.n.states      <- 2 # ***** NOTE ASSUMES Two states ONLY

  edge.mrf <- make.empty.field(adj.mat = edge.adj, parameterization.typ = "standard", plotQ = F)
  #edge.mrf <-train.mrf(edge.mrf, edge.samps.loc)

  # ***MLE parameters "by hand" instead of usinf train.mrf to get a little more control over the convergence.
  # Sometimes it takes a few runs of optim to get the gradient down to a reasonable size.

  # Auxiliary, gradient convenience function for optim.
  edge.gradient <- function(par, crf, ...) { crf$gradient }
  edge.mrf$par.stat <- mrf.stat(edge.mrf, instances = edge.samps.loc)

  convergedQ    <- 0
  loc.max.miter <- 10 # Max number of times to re-run optim.
  loc.miter     <- 1  # Iteration number for runing optim
  while(loc.miter <= loc.max.miter) {
    edge.opt.info  <- stats::optim(  # optimize parameters
      par          = edge.mrf$par,   # theta
      fn           = negloglik,      # objective function (CRFutil version, not CRF version)
      gr           = edge.gradient,  # grad of obj func
      crf          = edge.mrf,       # passed to fn/gr
      samples      = edge.samps.loc, # passed to fn/gr
      infer.method = infer.exact,    # passed to fn/gr
      update.crfQ  = TRUE,           # passed to fn/gr
      method       = "L-BFGS-B",
      control      = list(factr=10, trace = 1, REPORT=1))
    mag.grad <- sqrt(sum(edge.mrf$gradient^2)) # Magnitude, gradient of negloglik after optim has finished

    if(edge.opt.info$convergence == 0) {
      print("-------------------------------")
      print(paste("Macro iteration:", loc.miter, edge.opt.info$message))
      print(paste("|grad| =", mag.grad))
    } else {
      print("************* optim did not converge! *************")
      print(paste("Macro iteration:", loc.miter, edge.opt.info$message))
      print(paste("|grad L| =", mag.grad))
      warning("optim not converged.......")
    }

    # Check convergence of the gradient (macro convergence)
    if(mag.grad <= 1e-4){
      convergedQ <- 1
      print("Macro iterations converged! See above for gradient.")
      break()
    } else {
      print(paste("Macro iteration:", loc.miter))
      loc.miter <- loc.miter + 1
    }

    if(loc.miter > loc.max.miter) {
      print("Max number of macro iterations reached without meeting |grad| criteria. Terminating optimization process!")
    }

  }

  return(edge.mrf)

  # Compare fit beliefs to empirical probs. Is Bel(A|B) ~~ Bel(A) ~~ Pr(A) etc??

}


#' marginal edge beliefs
#'
#' Get edge beliefs Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1) from fit MRF
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bels <- function(edge.mrf.obj, node.names = NULL, printQ=FALSE){

  # if(ncol(edge.samples) != 2) {
  #   stop("Input one edge (two nodes) only!")
  # }

  if(is.null(node.names)) {
    loc.node.names <- c("X1", "X2")
  } else {
    loc.node.names <- node.names
  }

  infered.edge.bels <- make.gRbase.beliefs(
    inference.obj = infer.exact(edge.mrf.obj),
    node.names    = loc.node.names,
    edge.mat      = edge.mrf.obj$edges,
    state.nmes    = c("1","2"))

  bel.x1x2 <- infered.edge.bels$edge.beliefs[[1]]
  bel.x1   <- infered.edge.bels$node.beliefs[[1]]
  bel.x2   <- infered.edge.bels$node.beliefs[[2]]

  bel.x1gx2 <- ar_div(bel.x1x2, bel.x2)
  bel.x2gx1 <- ar_div(bel.x1x2, bel.x1)

  if(printQ==TRUE){

    print("---------")
    print("Bel(X1,X2)")
    print("---------")
    print(bel.x1x2)
    print("=======================")

    print("---------")
    print("Bel(X1)")
    print("---------")
    print(bel.x1)
    print("=======================")

    print("---------")
    print("Bel(X2)")
    print("---------")
    print(bel.x2)
    print("=======================")

    print("---------")
    print("Bel(X1|X2)")
    print("---------")
    print(bel.x1gx2)
    print("=======================")

    print("---------")
    print("Bel(X2|X1)")
    print("---------")
    print(bel.x2gx1)
    print("=======================")
  }

  edg.bel.info <- list(
    bel.x1x2,
    bel.x1,
    bel.x2,
    bel.x1gx2,
    bel.x2gx1
  )

  names(edg.bel.info) <- c(
    paste0("Bel(",loc.node.names[1],",",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[1],")"),
    paste0("Bel(",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[1],"|",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[2],"|",loc.node.names[1],")")
  )

  return(edg.bel.info)

}


#' marginal.edge.emp.pr
#'
#' Empirical Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1) estimates from edge data
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.emp.pr <- function(edge.samples, printQ=FALSE){

  if(ncol(edge.samples) != 2) {
    stop("Input one edge (two nodes) only!")
  }

  # Construct contingency table for edge data
  edge.samps.loc           <- edge.samples # Make a copy
  if(is.null(colnames(edge.samples))){
    colnames(edge.samps.loc) <- c("X1","X2")
  } else {
    colnames(edge.samps.loc) <- colnames(edge.samples)
  }
  X.edge.contingency       <- xtabs(~., data=edge.samps.loc)
  if(printQ==TRUE){
    print("-----------------------")
    print("Edge Contingency Table:")
    print("-----------------------")
    print(X.edge.contingency)
    print("=======================")
  }

  # Compute empirical marginal edge probabilities: Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1)
  nnmes    <- colnames(edge.samps.loc)
  pr.x1x2  <- X.edge.contingency/sum(X.edge.contingency)
  pr.x1    <- ar_marg(pr.x1x2, marg = nnmes[1])
  pr.x2    <- ar_marg(pr.x1x2, marg = nnmes[2])
  pr.x1gx2 <- ar_div(pr.x1x2, pr.x2)
  pr.x2gx1 <- ar_div(pr.x1x2, pr.x1)

  if(printQ==TRUE){

    print("---------")
    print("Pr(X1,X2)")
    print("---------")
    print(pr.x1x2)
    print("=======================")

    print("---------")
    print("Pr(X1)")
    print("---------")
    print(pr.x1)
    print("=======================")

    print("---------")
    print("Pr(X2)")
    print("---------")
    print(pr.x2)
    print("=======================")

    print("---------")
    print("Pr(X1|X2)")
    print("---------")
    print(pr.x1gx2)
    print("=======================")

    print("---------")
    print("Pr(X2|X1)")
    print("---------")
    print(pr.x2gx1)
    print("=======================")
  }

  edg.emp.pr.info <- list(
    X.edge.contingency,
    pr.x1x2,
    pr.x1,
    pr.x2,
    pr.x1gx2,
    pr.x2gx1
  )

  names(edg.emp.pr.info) <- c(
    "edge.contingency.tbl",
    paste0("Pr(",nnmes[1],",",nnmes[2],")"),
    paste0("Pr(",nnmes[1],")"),
    paste0("Pr(",nnmes[2],")"),
    paste0("Pr(",nnmes[1],"|",nnmes[2],")"),
    paste0("Pr(",nnmes[2],"|",nnmes[1],")")
  )

  return(edg.emp.pr.info)

}
