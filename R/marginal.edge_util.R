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
    edge.glm,               # fit glm object
    modl.coef.info[,1],     # log potential coeffs (energies) from glm fit
    glm.theta.est.vec,      # log potential coeffs (energies) according to p-values
    glm.node.pot,           # potentials info
    glm.node.pot.rescaled,
    glm.edge.pot,
    glm.edge.pot.rescaled
  )

  names(glm.pot.est.info) <- c(
    "glm.model",
    "glm.theta.raw",
    "glm.theta.est",
    "glm.poi.node.pot",
    "glm.poi.rescaled.node.pot",
    "glm.poi.edge.pot",
    "glm.poi.rescaled.edge.pot"
  )

  return(glm.pot.est.info)

}


#' marginal edge fit using vanalla baysian poisson regression with rstanarm
#'
#' XXXXXXX
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bayes.loglin <- function(edge.samples, prior.sd=NULL, prior_intercept.sd=NULL, stan.iter=NULL, prob.level=0.95, printQ=FALSE){

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

  # To not go too banannas, asume weakly informative normal priors
  # Set up width of normal priors
  if(is.null(prior.sd)) {
    prior.sd.loc <- 3
  } else {
    prior.sd.loc <- prior.sd
  }
  if(is.null(prior_intercept.sd)) {
    prior_intercept.sd.loc <- 3
  } else {
    prior_intercept.sd.loc <- prior.sd
  }

  if(is.null(stan.iter)) {
    stan.iter.loc <- 4000
  } else {
    stan.iter <- stan.iter.loc
  }

  # Use pre-coded vanalla poison regression in rstanarm:
  edge.bglm <-  stan_glm(Freq ~ X1 + X2 + omeg, # Don't worry. It will put in the intercept. Works like glm
                        family          = poisson,
                        data            = edge.data,
                        prior           = normal(0,prior.sd.loc),
                        prior_intercept = normal(0,prior_intercept.sd.loc),
                        chains          = 4,
                        iter            = stan.iter.loc,
                        cores           = 4)

  modl.summary <- summary(edge.bglm)
  if(printQ==TRUE){
    print(modl.summary)
  }

  modl.coef.info  <- edge.bglm$coefficients # parameter medians
  post.param.samp <- as.matrix(edge.bglm)
  alp  <- post.param.samp[,1]               # posterior parameter samples
  tau1 <- post.param.samp[,2]
  tau2 <- post.param.samp[,3]
  omeg <- post.param.samp[,4]

  alp.int    <- HPDI(samples = alp,  prob = prob.level) # parameter prob intervals
  tau1.int   <- HPDI(samples = tau1, prob = prob.level)
  tau2.int   <- HPDI(samples = tau2, prob = prob.level)
  omeg.int   <- HPDI(samples = omeg, prob = prob.level)
  param.ints <- rbind(alp.int, tau1.int, tau2.int, omeg.int)

  modl.coef.info           <- data.frame(modl.coef.info,
                                         param.ints,
                                         !cbind( param.ints[,1] <= 0 & param.ints[,2] >= 0 ) )
  bglm.par.est             <- modl.coef.info[,4]*modl.coef.info[,1]
  modl.coef.info           <- data.frame(modl.coef.info, bglm.par.est)

  colnames(modl.coef.info) <- c("param.median", names(alp.int), "0.NOT.covered?", "bglm.pars")
  rownames(modl.coef.info) <- c("(Intercept)", "tau1", "tau2", "omega")
  print(modl.coef.info)

  # Transform and re-scale posterior sample of parameters to sample of potentials
  pot.tau1 <- exp(tau1)/exp(-tau1)
  pot.tau2 <- exp(tau2)/exp(-tau2)
  pot.omeg <- exp(omeg)/exp(-omeg)

  raw.theta                         <- modl.coef.info[,1]
  names(raw.theta)                  <- c("(Intercept)", "tau1", "tau2", "omega")
  med.theta                         <- bglm.par.est[-1]
  names(med.theta)                  <- c("tau1", "tau2", "omega")
  rescaled.posterior.pots           <- cbind(pot.tau1, pot.tau2, pot.omeg)
  colnames(rescaled.posterior.pots) <- c("pot.tau1","pot.tau2","pot.omega")

  edge.bglm.info <- list(
    edge.bglm,              # rstanarm fit object. Contains unscaled original posterior samples
    raw.theta,              # all log coefs (energies) from stan fit
    med.theta,              # essential and cleaned up log coefs (energies) according to HPDI
    rescaled.posterior.pots # potentials info
  )

  names(edge.bglm.info) <- c(
    "rstanarm.obj",
    "raw.theta",
    "theta.par.meds",
    "rescaled.posterior.pots"
  )

  return(edge.bglm.info)

}


#' marginal edge fit with MLE within the CRF framework
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

}


#' Marginal edge beliefs from parameter point estimates.
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


#' Marginal edge beliefs from bayes posterior sample
#'
#' Get edge beliefs Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1) from bayes loglin fit
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bels.bayes <- function(post.edge.param.samp, node.names = NULL, state.names = NULL, printQ=FALSE){

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

  loc.graph.eq   <- ~X1loc + X2loc + X1loc:X2loc # Marginal edge graph: X1--X2
  edge.mrf.obj   <- make.empty.field(graph.eq = loc.graph.eq, parameterization.typ = "standard", plotQ = F)
  post.bel.x1x2  <- array(0,c(nrow(post.edge.param.samp), 4))
  post.bel.x1    <- array(0,c(nrow(post.edge.param.samp), 2))
  post.bel.x2    <- array(0,c(nrow(post.edge.param.samp), 2))
  post.bel.x1gx2 <- array(0,c(nrow(post.edge.param.samp), 4))
  post.bel.x2gx1 <- array(0,c(nrow(post.edge.param.samp), 4))

  # For each posterior potential sample compute Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1)
  for(i in 1:nrow(post.edge.param.samp)) {

    post.pots <- post.edge.param.samp[i,]

    #TMP
    # print(paste("================== SAMPLE:",i,"======================"))
    # print("Posterior potential vector:")
    # print(post.pots)
    # print("---------")


    edge.mrf.obj$node.pot <- rbind(
      c(post.pots[1], 1),
      c(post.pots[2], 1)
    )
    edge.mrf.obj$edge.pot[[1]] <- rbind(
      c(post.pots[3], 1           ),
      c(1           , post.pots[3])
    )

    #TMP
    # print("Input node potential matrix:")
    # print(edge.mrf.obj$node.pot)
    # print("---------")
    #
    # #TMP
    # print("Input edge potential matrix:")
    # print(edge.mrf.obj$edge.pot)
    # print("---------")

    # Get the edge beliefs and flatten them into stackable format    PROBLEM HERE?????????????
    infered.edge.bels <- marginal.edge.bels(
      edge.mrf.obj = edge.mrf.obj,
      node.names   = loc.node.names,
      state.names  = loc.state.names)

    #TMP
    # print("Infered edge beliefs:")
    # print(infered.edge.bels)
    # print("---------")

    infered.edge.bels <- flatten.marginal.edge.beliefs(
      marginal.edge.beliefs.list = infered.edge.bels)

    #TMP
    # print("Flattened Infered edge beliefs:")
    # print(infered.edge.bels)
    # print("---------")

    bel.x1x2  <- infered.edge.bels[[1]]
    bel.x1    <- infered.edge.bels[[2]]
    bel.x2    <- infered.edge.bels[[3]]
    bel.x1gx2 <- infered.edge.bels[[4]]
    bel.x2gx1 <- infered.edge.bels[[5]]

    post.bel.x1x2[i,]  <- bel.x1x2
    post.bel.x1[i,]    <- bel.x1
    post.bel.x2[i,]    <- bel.x2
    post.bel.x1gx2[i,] <- bel.x1gx2
    post.bel.x2gx1[i,] <- bel.x2gx1

    #print(i)
  }

  # Using the last iteration info get the column names for each sample matrix
  colnames(post.bel.x1x2)  <- colnames(infered.edge.bels[[1]])
  colnames(post.bel.x1)    <- colnames(infered.edge.bels[[2]])
  colnames(post.bel.x2)    <- colnames(infered.edge.bels[[3]])
  colnames(post.bel.x1gx2) <- colnames(infered.edge.bels[[4]])
  colnames(post.bel.x2gx1) <- colnames(infered.edge.bels[[5]])

  post.edg.bel.info <- list(
    post.bel.x1x2,
    post.bel.x1,
    post.bel.x2,
    post.bel.x1gx2,
    post.bel.x2gx1
  )

  names(post.edg.bel.info) <- c(
    paste0("Bel(",loc.node.names[1],",",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[1],")"),
    paste0("Bel(",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[1],"|",loc.node.names[2],")"),
    paste0("Bel(",loc.node.names[2],"|",loc.node.names[1],")")
  )

  return(post.edg.bel.info)

}


#' Marginal edge empirical probs from a sample of data
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
    print("Pr(X1|X2)")
    print("---------")
    print(pr.x1gx2)

    print("=======================")
    print("---------")
    print("Pr(X1)")
    print("---------")
    print(pr.x1)
    print("=======================")

    print("---------")
    print("Pr(X2|X1)")
    print("---------")
    print(pr.x2gx1)
    print("=======================")

    print("---------")
    print("Pr(X2)")
    print("---------")
    print(pr.x2)
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


#' Flatten belief matrices output by marginal.edge.bels to make them stackable
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
flatten.marginal.edge.beliefs <- function(marginal.edge.beliefs.list, printQ=FALSE){

  bel.nmes <- names(marginal.edge.beliefs.list)
  #print(bel.nmes)

  # Flatten Bel(X1,X2)
  tij       <- marginal.edge.beliefs.list[[ bel.nmes[1] ]]              # mij Belief table
  stij      <- colnames(tij)                                            # state names
  din       <- dimnames( marginal.edge.beliefs.list[[ bel.nmes[1] ]] )  # dnij Node names for table
  nis       <- names(din)                                               # nij Node names collapsed
  x1x2.lbls <- array("",c(1,4))
  x1x2.vals <- array(NA,c(1,4))
  count <- 1
  for(i in 1:2){  # col
    for(j in 1:2){ # row
      x1x2.lbls[count] <- paste0(nis[2], "=", stij[i], ",",nis[1],"=",stij[j])
      x1x2.vals[count] <- tij[j,i]
      count <- count + 1
    }
  }
  colnames(x1x2.vals) <- x1x2.lbls
  if(printQ==TRUE) {
    print(bel.nmes[1])
    print(x1x2.vals)
  }

  #  Bel(X1)
  tij     <- marginal.edge.beliefs.list[[ bel.nmes[2] ]]              # mij Belief table
  stij    <- names(tij)                                               # state names
  din     <- dimnames( marginal.edge.beliefs.list[[ bel.nmes[2] ]] )  # dnij Node names for table
  nis     <- names(din)                                               # nij Node names collapsed
  x1.lbls <- array("",c(1,2))
  x1.vals <- array(NA,c(1,2))
  count <- 1
  for(i in 1:2){  # col
    x1.lbls[count] <- paste0(nis[1], "=", stij[i])
    x1.vals[count] <- tij[i]
    count <- count + 1
  }
  colnames(x1.vals) <- x1.lbls
  if(printQ==TRUE) {
    print(bel.nmes[2])
    print(x1.vals)
  }

  #  Bel(X2)
  tij     <- marginal.edge.beliefs.list[[ bel.nmes[3] ]]              # mij Belief table
  stij    <- names(tij)                                               # state names
  din     <- dimnames( marginal.edge.beliefs.list[[ bel.nmes[3] ]] )  # dnij Node names for table
  nis     <- names(din)                                               # nij Node names collapsed
  x2.lbls <- array("",c(1,2))
  x2.vals <- array(NA,c(1,2))
  count <- 1
  for(i in 1:2){  # col
    x2.lbls[count] <- paste0(nis[1], "=", stij[i])
    x2.vals[count] <- tij[i]
    count <- count + 1
  }
  colnames(x2.vals) <- x2.lbls
  if(printQ==TRUE) {
    print(bel.nmes[3])
    print(x2.vals)
  }

  # Flatten Bel(X1|X2)
  tij        <- marginal.edge.beliefs.list[[ bel.nmes[4] ]]              # mij Belief table
  stij       <- colnames(tij)                                            # state names
  din        <- dimnames( marginal.edge.beliefs.list[[ bel.nmes[4] ]] )  # dnij Node names for table
  nis        <- names(din)                                               # nij Node names collapsed
  x1gx2.lbls <- array("",c(1,4))
  x1gx2.vals <- array(NA,c(1,4))
  count <- 1
  for(i in 1:2){  # col
    for(j in 1:2){ # row
      x1gx2.lbls[count] <- paste0(nis[2], "=", stij[i], "|", nis[1],"=", stij[j])
      x1gx2.vals[count] <- tij[j,i]
      count <- count + 1
    }
  }
  colnames(x1gx2.vals) <- x1gx2.lbls
  if(printQ==TRUE) {
    print(bel.nmes[4])
    print(x1gx2.vals)
  }

  # Flatten Bel(X2|X1)
  tij        <- marginal.edge.beliefs.list[[ bel.nmes[5] ]]              # mij Belief table
  stij       <- colnames(tij)                                            # state names
  din        <- dimnames( marginal.edge.beliefs.list[[ bel.nmes[5] ]] )  # dnij Node names for table
  nis        <- names(din)                                               # nij Node names collapsed
  x2gx1.lbls <- array("",c(1,4))
  x2gx1.vals <- array(NA,c(1,4))
  count <- 1
  for(i in 1:2){  # col
    for(j in 1:2){ # row
      x2gx1.lbls[count] <- paste0(nis[2], "=", stij[i], "|", nis[1], "=", stij[j])
      x2gx1.vals[count] <- tij[j,i]
      count <- count + 1
    }
  }
  colnames(x2gx1.vals) <- x2gx1.lbls
  if(printQ==TRUE) {
    print(bel.nmes[5])
    print(x2gx1.vals)
  }

  flat.bels <- list(
    x1x2.vals,
    x1.vals,
    x2.vals,
    x1gx2.vals,
    x2gx1.vals
  )

  names(flat.bels) <- c(
    bel.nmes[1],
    bel.nmes[2],
    bel.nmes[3],
    bel.nmes[4],
    bel.nmes[5]
  )

  return(flat.bels)

}
