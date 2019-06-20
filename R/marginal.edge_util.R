#' marginal edge fit using poisson regression glm
#'
#' XXXXXXX
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.loglin <- function(edge.samples, conf.level=0.95, printQ=FALSE){

  if(ncol(edge.samples) != 2) {
    stop("Input one edge (two nodes) only!")
  }

  # Construct contingency table for edge data. For now just use X1 and X2 as node names
  edge.samps.loc           <- edge.samples # Make a copy
  # if(is.null(colnames(edge.samples))){
  #   colnames(edge.samps.loc) <- c("X1","X2")
  # } else {
  #   colnames(edge.samps.loc) <- colnames(edge.samples)
  # }
  colnames(edge.samps.loc) <- c("X1","X2")
  edge.contingency         <- xtabs(~., data=edge.samps.loc)
  edge.freq.table          <- as.data.frame(edge.contingency)
  edge.freq                <- edge.freq.table[,3]

  # Fit a two node loglinear (Poisson) model to examine the edge potentials and resulting node/edge beliefs.

  #Contrasts sum model matrix for X1--X2
  modl.mat <- rbind(
    c(1,  1,  1,  1),
    c(1, -1,  1, -1),
    c(1,  1, -1, -1),
    c(1, -1, -1,  1))
  colnames(modl.mat)       <- c("(Intercept)", "X1", "X2", "X1:X2")
  edge.glm                 <- glm(edge.freq ~ modl.mat[,1] + modl.mat[,2] + modl.mat[,3] + modl.mat[,4] -1, family = poisson(link="log"))
  modl.summary             <- summary(edge.glm)
  if(printQ==TRUE){
    print(modl.summary)
  }

  modl.coef.info           <- modl.summary$coefficients
  modl.coef.info           <- data.frame(modl.coef.info, modl.coef.info[,4] < (1-conf.level))
  modl.coef.info           <- modl.coef.info[,-3]
  glm.par.est              <- modl.coef.info[,4]*modl.coef.info[,1]
  modl.coef.info           <- data.frame(modl.coef.info,glm.par.est)

  colnames(modl.coef.info) <- c("theta.hat","std.err","p.val","Reject.H0:theta=0?","glm.pars")
  rownames(modl.coef.info) <- colnames(modl.mat)
  print(modl.coef.info)

  # Put theta est into potential matrix format and re-scale:
  glm.theta.est.vec <- glm.par.est[-1]

  glm.node.pot <- rbind(
    c( exp(glm.theta.est.vec[1]), exp(-glm.theta.est.vec[1]) ),
    c( exp(glm.theta.est.vec[2]), exp(-glm.theta.est.vec[2]) )
  )

  glm.edge.pot <- rbind(
    c(exp(glm.theta.est.vec[3]),  exp(-glm.theta.est.vec[3])),
    c(exp(-glm.theta.est.vec[3]), exp(glm.theta.est.vec[3]))
  )

  # Re-scale node pot matrix wrt second column elements
  glm.node.pot.rescaled     <- glm.node.pot
  glm.node.pot.rescaled[1,] <- glm.node.pot.rescaled[1,]/glm.node.pot.rescaled[1,2]
  glm.node.pot.rescaled[2,] <- glm.node.pot.rescaled[2,]/glm.node.pot.rescaled[2,2]

  # Re-scale edge pot matrix wrt off-diagonal elements
  glm.edge.pot.rescaled     <- glm.edge.pot
  glm.edge.pot.rescaled[1,] <- glm.edge.pot.rescaled[1,]/glm.edge.pot.rescaled[1,2]
  glm.edge.pot.rescaled[2,] <- glm.edge.pot.rescaled[2,]/glm.edge.pot.rescaled[2,1]

  glm.pot.est.info <- list(
    modl.coef.info[,1],     # log potential coeffs (energies) from glm fit
    glm.theta.est.vec,      # log potential coeffs (energies) according to p-values
    glm.node.pot,           # potentials info
    glm.node.pot.rescaled,
    glm.edge.pot,
    glm.edge.pot.rescaled
  )

  names(glm.pot.est.info) <- c(
    "glm.theta.raw",
    "glm.theta.est",
    "glm.poi.node.pot",
    "glm.poi.rescaled.node.pot",
    "glm.poi.edge.pot",
    "glm.poi.rescaled.edge.pot"
  )

  return(glm.pot.est.info)

}


#' marginal edge fit using poisson regression glm
#'
#' XXXXXXX
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bayes.loglin <- function(edge.samples, printQ=FALSE){

  if(ncol(edge.samples) != 2) {
    stop("Input one edge (two nodes) only!")
  }

  # Construct contingency table for edge data. For now just use X1 and X2 as node names
  edge.samps.loc           <- edge.samples # Make a copy
  colnames(edge.samps.loc) <- c("X1","X2")
  edge.contingency         <- xtabs(~., data=edge.samps.loc)
  edge.freq.table          <- as.data.frame(edge.contingency)
  edge.freq                <- edge.freq.table[,3]

  # Fit a two node loglinear (Poisson) model to examine the edge potentials and resulting node/edge beliefs.

  #Contrasts sum model matrix for X1--X2
  modl.mat <- rbind(
    c(1,  1,  1,  1),
    c(1, -1,  1, -1),
    c(1,  1, -1, -1),
    c(1, -1, -1,  1))
  edge.data <- data.frame(modl.mat, edge.freq)
  colnames(edge.data)        <- c("alp", "X1", "X2", "omeg","Freq") # omeg is the edge param, alp is the intercept
  #print(edge.data)
  edge.glm <-  stan_glm(Freq ~ X1 + X2 + omeg,
                        family = poisson,
                        data = edge.data,
                        prior = normal(0,3),
                        prior_intercept = normal(0,10),
                        chains = 4,
                        cores = 4)
  return(edge.glm)
  # modl.summary             <- summary(edge.glm)
  # if(printQ==TRUE){
  #   print(modl.summary)
  # }
  #
  # modl.coef.info           <- modl.summary$coefficients
  # modl.coef.info           <- data.frame(modl.coef.info, modl.coef.info[,4] < (1-conf.level))
  # modl.coef.info           <- modl.coef.info[,-3]
  # glm.par.est              <- modl.coef.info[,4]*modl.coef.info[,1]
  # modl.coef.info           <- data.frame(modl.coef.info,glm.par.est)
  #
  # colnames(modl.coef.info) <- c("theta.hat","std.err","p.val","Reject.H0:theta=0?","glm.pars")
  # rownames(modl.coef.info) <- colnames(modl.mat)
  # print(modl.coef.info)
  #
  # # Put theta est into potential matrix format and re-scale:
  # glm.theta.est.vec <- glm.par.est[-1]
  #
  # glm.node.pot <- rbind(
  #   c( exp(glm.theta.est.vec[1]), exp(-glm.theta.est.vec[1]) ),
  #   c( exp(glm.theta.est.vec[2]), exp(-glm.theta.est.vec[2]) )
  # )
  #
  # glm.edge.pot <- rbind(
  #   c(exp(glm.theta.est.vec[3]),  exp(-glm.theta.est.vec[3])),
  #   c(exp(-glm.theta.est.vec[3]), exp(glm.theta.est.vec[3]))
  # )
  #
  # # Re-scale node pot matrix wrt second column elements
  # glm.node.pot.rescaled     <- glm.node.pot
  # glm.node.pot.rescaled[1,] <- glm.node.pot.rescaled[1,]/glm.node.pot.rescaled[1,2]
  # glm.node.pot.rescaled[2,] <- glm.node.pot.rescaled[2,]/glm.node.pot.rescaled[2,2]
  #
  # # Re-scale edge pot matrix wrt off-diagonal elements
  # glm.edge.pot.rescaled     <- glm.edge.pot
  # glm.edge.pot.rescaled[1,] <- glm.edge.pot.rescaled[1,]/glm.edge.pot.rescaled[1,2]
  # glm.edge.pot.rescaled[2,] <- glm.edge.pot.rescaled[2,]/glm.edge.pot.rescaled[2,1]
  #
  # glm.pot.est.info <- list(
  #   modl.coef.info[,1],     # log potential coeffs (energies) from glm fit
  #   glm.theta.est.vec,      # log potential coeffs (energies) according to p-values
  #   glm.node.pot,           # potentials info
  #   glm.node.pot.rescaled,
  #   glm.edge.pot,
  #   glm.edge.pot.rescaled
  # )
  #
  # names(glm.pot.est.info) <- c(
  #   "glm.theta.raw",
  #   "glm.theta.est",
  #   "glm.poi.node.pot",
  #   "glm.poi.rescaled.node.pot",
  #   "glm.poi.edge.pot",
  #   "glm.poi.rescaled.edge.pot"
  # )
  #
  # return(glm.pot.est.info)

}


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

  edge.samps.loc <- edge.samples # Make a copy
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
marginal.edge.bels <- function(edge.mrf.obj, node.names = NULL, state.names = NULL, printQ=FALSE){

  # if(ncol(edge.samples) != 2) {
  #   stop("Input one edge (two nodes) only!")
  # }

  if(is.null(node.names)) {
    loc.node.names <- c("X1", "X2")
  } else {
    loc.node.names <- node.names
  }

  if(is.null(state.names)) {
    loc.state.names <- c("1","2")
  } else {
    loc.state.names <- state.names
  }

  infered.edge.bels <- make.gRbase.beliefs(
    inference.obj = infer.exact(edge.mrf.obj),
    node.names    = loc.node.names,
    edge.mat      = edge.mrf.obj$edges,
    state.nmes    = loc.state.names)

  bel.x1x2 <- infered.edge.bels$edge.beliefs[[1]]
  bel.x1   <- infered.edge.bels$node.beliefs[[1]]
  bel.x2   <- infered.edge.bels$node.beliefs[[2]]

  bel.x1gx2 <- ar_div(bel.x1x2, bel.x2)
  bel.x2gx1 <- ar_div(bel.x1x2, bel.x1)

  if(printQ==TRUE){

    print("----------")
    print("Bel(X1,X2)")
    print("----------")
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

    print("----------")
    print("Bel(X1|X2)")
    print("----------")
    print(bel.x1gx2)
    print("=======================")

    print("----------")
    print("Bel(X2|X1)")
    print("----------")
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


#' marginal edge empirical probs
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
