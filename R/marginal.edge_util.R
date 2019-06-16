#' marginal.edge.info
#'
#' XXXXXXX
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.info <- function(edge.samples){

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
  print(X.edge.contingency)

  # Compute empirical marginal edge probabilities: Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1)
  nnmes    <- colnames(edge.samps.loc)
  pr.x1x2  <- X.edge.contingency/sum(X.edge.contingency)
  pr.x1    <- ar_marg(pr.x1x2, marg = nnmes[1])
  pr.x2    <- ar_marg(pr.x1x2, marg = nnmes[2])
  pr.x1gx2 <- ar_div(pr.x1x2, pr.x2)
  pr.x2gx1 <- ar_div(pr.x1x2, pr.x1)

  print("Pr(X1,X2)")
  print(pr.x1x2)
  print("Pr(X1)")
  print(pr.x1)
  print("Pr(X2)")
  print(pr.x2)
  print("Pr(X1|X2)")
  print(pr.x1gx2)
  print("Pr(X2|X1)")
  print(pr.x2gx1)

  # Fit a two node MRF to examine the edge potentials and resulting node/edge beliefs.
  edge.grphf    <- ~EX1:EX2
  edge.adj      <- ug(edge.grphf, result="matrix")
  colnames(edge.adj) <- nnmes
  rownames(edge.adj) <- nnmes
  edge.n.nodes  <- 2
  edge.n.states <- 2 # ***** NOTE ASSUMES Two states ONLY

  edge.mrf <- make.empty.field(adj.mat = edge.adj, parameterization.typ = "standard", plotQ = F)
  #edge.mrf <-train.mrf(edge.mrf, edge.samps.loc)

  # Auxiliary, gradient convenience function for optim.
  gradient <- function(par, crf, ...) { crf$gradient }
  edge.mrf$par.stat <- mrf.stat(edge.mrf, instances = edge.samps.loc)

  edge.opt.info  <- stats::optim(  # optimize parameters
    par          = edge.mrf$par,   # theta
    fn           = negloglik,      # objective function
    gr           = gradient,       # grad of obj func
    crf          = edge.mrf,       # passed to fn/gr
    samples      = edge.samps.loc, # passed to fn/gr
    infer.method = infer.exact,    # passed to fn/gr
    update.crfQ  = TRUE,           # passed to fn/gr
    method       = "L-BFGS-B",
    control      = list(factr=10, trace = 1, REPORT=1))
  mag.grad <- sqrt(sum(edge.mrf$gradient^2)) # Magnitude, gradient of negloglik after optim has finished

  if(edge.opt.info$convergence == 0) {
    print(paste("Macro iteration:", edge.opt.info$message))
    print(paste("|grad| =", mag.grad))
  } else {
    warning("optim not converged.......")
    print("Macro iteration:", edge.opt.info$message)
    print("|grad L| =", mag.grad)
  }

  #print(edge.mrf$par)
  #print(edge.mrf$gradient)
  return(edge.mrf)



  # Compare fit beliefs to empirical probs. Is Bel(A|B) ~~ Bel(A) ~~ Pr(A) etc??

}
