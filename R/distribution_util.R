#' Compute the joint distribution from the node and edge energies
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
distribution.from.energies <- function(state.space, edges.mat, node.energies, edge.energies, energy.func, ff){

  num.states      <- nrow(state.space)
  state.energies  <- sapply(1:num.states, function(xx){energy.func(state.space[xx,], edges.mat, node.energies, edge.energies, ff)})
  logZZ           <- logsumexp2(state.energies)
  log.state.probs <- state.energies-logZZ

  dist.info <- list(exp(log.state.probs), logZZ)
  names(dist.info) <- c("state.probs", "logZ")
  return(dist.info)

}


#' Compute the joint distribution from the node and edge potentials using gRbase functions
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will state probalities in gRbase table form as well as logZ
#'
#'
#' @export
distribution.from.potentials <- function(gRbase.node.potentials, gRbase.edge.potentials){

  num.nodes <- length(gRbase.node.potentials)
  num.edges <- length(gRbase.edge.potentials)

  prod.node.pots <- tableMult(gRbase.node.potentials[[2]], gRbase.node.potentials[[1]])
  if(num.nodes > 2){
    for(i in 3:num.nodes){
      prod.node.pots <- tableMult(prod.node.pots, gRbase.node.potentials[[i]])
    }
  }

  #prod.edge.pots <- tableMult(gRbase.edge.potentials[[2]], gRbase.edge.potentials[[1]])
  if(num.edges > 2){
    prod.edge.pots <- tableMult(gRbase.edge.potentials[[2]], gRbase.edge.potentials[[1]])
    for(i in 3:num.edges){
      prod.edge.pots <- tableMult(prod.edge.pots, gRbase.edge.potentials[[i]])
    }
  } else { # For only two nodes there is only one edge
    prod.edge.pots <- gRbase.edge.potentials[[1]]
  }

  # Direct normalization:
  #state.probs <- tableMult(prod.edge.pots, prod.node.pots)
  #ZZ <- sum(state.probs)
  #state.probs <- state.probs/ZZ

  # Assume the prod pots can get a little rowdy. Normalize on log scale instead:
  log.state.prod.pots <- log(tableMult(prod.edge.pots, prod.node.pots))
  logZZ               <- logsumexp2(log.state.prod.pots)
  log.state.probs     <- log.state.prod.pots - logZZ

  dist.info <- list(exp(log.state.probs), logZZ)
  names(dist.info) <- c("state.probs", "logZ")
  return(dist.info)

}

#' Compute the joint distribution, assuming a crf object is passed in with gRbase formatted potentials and/or energies
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will state probalities in gRbase table form as well as logZ
#'
#'
#' @export
distribution.calc <- function(crf, logZ.calc=NULL){  # change logZ.calc to inference.meth at some point

  num.nodes <- crf$n.nodes
  num.edges <- crf$n.edges

  prod.node.pots <- tableMult(crf$node.potentials[[2]], crf$node.potentials[[1]])
  if(num.nodes > 2){
    for(i in 3:num.nodes){
      prod.node.pots <- tableMult(prod.node.pots, crf$node.potentials[[i]])
    }
  }

  #prod.edge.pots <- tableMult(gRbase.edge.potentials[[2]], gRbase.edge.potentials[[1]])
  if(num.edges > 2){
    prod.edge.pots <- tableMult(crf$edge.potentials[[2]], crf$edge.potentials[[1]])
    for(i in 3:num.edges){
      prod.edge.pots <- tableMult(prod.edge.pots, crf$edge.potentials[[i]])
    }
  } else { # For only two nodes there is only one edge
    prod.edge.pots <- crf$edge.potentials[[1]]
  }

  # Direct normalization:
  #state.probs <- tableMult(prod.edge.pots, prod.node.pots)
  #ZZ <- sum(state.probs)
  #state.probs <- state.probs/ZZ

  # Assume the prod pots can get a little rowdy. Normalize on log scale instead:
  log.state.prod.pots <- log(tableMult(prod.edge.pots, prod.node.pots)) # log un-normalized joint dist.

  if(is.null(logZ.calc)){
    logZZ <- logsumexp2(log.state.prod.pots)
  } else {
    logZZ <- logZ.calc(crf)$logZ # Use one of the inference methods implemented in CRF instead
  }

  log.state.probs     <- log.state.prod.pots - logZZ

  dist.info <- list(exp(log.state.probs), logZZ)
  names(dist.info) <- c("state.probs", "logZ")
  return(dist.info)

}


#' Compute the joint distribution, with just an input crf object
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#'
#' @details A little shorter version and more generic version of distribution.calc above.
#'
#' @return The function will compute all configuration probabilities as well as logZ
#'
#'
#' @export
compute.full.distribution <- function(crf.obj, node.names=NULL, state.names=NULL){

  if(is.null(node.names)) {
    node.names.loc <- as.character(1:crf.obj$n.nodes)
  }

  if(is.null(state.names)) {
    state.names.loc <- c("1", "2")
  }

  pot.info.loc     <- make.gRbase.potentials(crf.obj, node.names = node.names.loc, state.nmes = state.names.loc)
  gR.dist.info.loc <- distribution.from.potentials(pot.info.loc$node.potentials, pot.info.loc$edge.potentials)

  logZ.loc             <- gR.dist.info.loc$logZ
  joint.dist.info.loc  <- as.data.frame(as.table(gR.dist.info.loc$state.probs))
  joint.dist.info.list <- list(joint.dist.info.loc, logZ.loc)
  names(joint.dist.info.list) <- c("joint.distribution", "logZ")

  return(joint.dist.info.list)

}


#' Align joint distributions in table form, to a reference joint distribution
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#'
#' @details All joint distributions need to have the same numbers of nodes and configurations
#'
#' @return The function will compute all configuration probabilities as well as logZ
#'
#'
#' @export
align.distributions <- function(ref.dist.frame, dist.frame.list, include.configsQ=T){

  # Assumes last column is the distribution
  num.nodes.ref  <- ncol(ref.dist.frame) - 1
  node.names.ref <- colnames(ref.dist.frame)[1:num.nodes.ref]

  # Number of reference confgs:
  num.configs.ref <- nrow(ref.dist.frame)

  dist.probs <- array(-1.0, c(nrow(ref.dist.frame), length(dist.frame.list)))
  for(i in 1:length(dist.frame.list)) {

    dist.i             <- dist.frame.list[[i]]
    num.nodes.dist.i   <- ncol(dist.i) - 1
    node.names.dist.i  <- colnames(dist.i)[1:num.nodes.dist.i]
    num.configs.dist.i <- nrow(dist.i)

    # Make sure number of nodes are the same between dists
    if(num.nodes.dist.i != num.nodes.ref) {
      stop(paste("Dist ", i, " and Ref do not have the same number of nodes!"))
    }

    # Make sure number of configurations are the same between dists
    if(num.configs.dist.i != num.configs.ref) {
      stop(paste("Dist ", i, " and Ref do not have the same number of configurations!"))
    }

    # Permute the columns of dist.i to match the reference dist.:
    col.perm.i <- sapply(1:num.nodes.ref, function(xx){which(colnames(dist.i) == node.names.ref[xx])})
    dist.i     <- dist.i[, c(col.perm.i, ncol(dist.i))]

    # Rearrange state indices to be in the same order between the two models:
    rearr.idxs <- sapply(1:nrow(ref.dist.frame),function(xx){row.match(ref.dist.frame[xx,1:num.nodes.ref], table = dist.i[,1:num.nodes.dist.i])})
    dist.i     <- dist.i[rearr.idxs,]
    #print(data.frame(ref.dist.frame, dist.i))

    dist.probs[,i] <- dist.i[, ncol(dist.i)]

  }

  # Tack on configurations if requested
  if(include.configsQ == T) {
    dist.probs <- data.frame(dist.i[, 1:(ncol(ref.dist.frame)-1) ], dist.probs)
  }

  return(dist.probs)

}


#' Compute the pseudolikelihoods from the node and edge energies. Assumes only 2 states/node
#' Pr(X) ~= Prod Pr(Xi | X/Xi) Besag 1975
#'
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
pseudolikelihoods.from.energies <- function(state.space, adjacent.nodes, edges.mat, node.energies, edge.energies, conditional.energy.func, ff){


  num.nodes       <- ncol(state.space)
  num.states      <- nrow(state.space)

  condtional.energies            <- array(NA, c(num.states, num.nodes))
  complement.condtional.energies <- array(NA, c(num.states, num.nodes))
  conditional.Prs                <- array(NA, c(num.states, num.nodes))
  complement.conditional.Prs     <- array(NA, c(num.states, num.nodes))
  pseudo.likelihoods             <- array(NA, num.states)
  conditional.Zs                 <- array(NA, c(num.states, num.nodes))
  for(i in 1:num.states) {

    # Conditional node energy E(Xi | X/Xi)
    node.condtional.energies <- sapply(1:num.nodes,function(xx){
      conditional.config.energy(state.space[i,],
                                condition.element.number=xx,
                                adj.node.list = adjacent.nodes,
                                edge.mat = edges.mat,
                                one.lgp = node.energies,
                                two.lgp = edge.energies,
                                ff = ff)})

    condtional.energies[i,] <- node.condtional.energies

    #  Complenent conditional node energy E(not-Xi | X/Xi)
    node.complement.condtional.energies <- sapply(1:num.nodes,function(xx){
      conditional.config.energy(complement.at.idx(state.space[i,],xx),
                                condition.element.number=xx,
                                adj.node.list = adjacent.nodes,
                                edge.mat = edges.mat,
                                one.lgp = node.energies,
                                two.lgp = edge.energies,
                                ff = ff)})

    complement.condtional.energies[i,] <- node.complement.condtional.energies

    # Zs for each conditional:
    node.conditional.Zs <- exp(node.condtional.energies) + exp(node.complement.condtional.energies)
    conditional.Zs[i,]  <- node.conditional.Zs

    # Pr(Xi | X/Xi)
    Prs.ce              <- exp(node.condtional.energies)/node.conditional.Zs
    conditional.Prs[i,] <- Prs.ce

    # Pr(not-Xi | X/Xi)
    #Prs.cce <- exp(node.complement.condtional.energies)/(exp(node.condtional.energies) + exp(node.complement.condtional.energies))
    Prs.cce                        <- 1 - Prs.ce
    complement.conditional.Prs[i,] <- Prs.cce

    # pseudo-liklihood = Prod Pr(Xi | X/Xi)
    pseudo.lik            <- prod(Prs.ce)
    pseudo.likelihoods[i] <- pseudo.lik

  }

  # Product pseudo likelihoods are not normalized wrt to \Pr({\bf X}), so lets just renormalize them
  # so atleast they sum to 1 (ie L1 normalize the product pseudo likelihoods)
  renormed.pseudo.likelihoods <- pseudo.likelihoods/sum(pseudo.likelihoods)

  dist.info <- list(
    condtional.energies,
    complement.condtional.energies,
    conditional.Zs,
    conditional.Prs,
    complement.conditional.Prs,
    pseudo.likelihoods,
    renormed.pseudo.likelihoods
  )
  names(dist.info) <- c(
    "condtional.energies",
    "complement.condtional.energies",
    "conditional.Zs",
    "conditional.Prs",
    "complement.conditional.Prs",
    "pseudo.likelihoods",
    "L1.renormalized.pseudo.likelihoods"
  )

  return(dist.info)

}


#' Revised: Compute the pseudolikelihoods from the node and edge energies. Assumes only 2 states/node
#' Pr(X) ~= Prod Pr(Xi | X/Xi) Besag 1975
#' The function will XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
pseudolikelihoods.from.energies2 <- function(state.space, crf, cond.en.form="feature", ff){

  num.nodes      <- ncol(state.space)
  num.states     <- nrow(state.space)
  #adjacent.nodes <- crf$adj.nodes
  #edges.mat      <- crf$edges
  theta.par      <- crf$par

  condtional.energies            <- array(NA, c(num.states, num.nodes))
  complement.condtional.energies <- array(NA, c(num.states, num.nodes))
  conditional.Prs                <- array(NA, c(num.states, num.nodes))
  complement.conditional.Prs     <- array(NA, c(num.states, num.nodes))
  pseudo.likelihoods             <- array(NA, num.states)
  conditional.Zs                 <- array(NA, c(num.states, num.nodes))

  if(cond.en.form=="feature.function"){
    cond.en.func <- conditional.config.energy
  } else {
    if(cond.en.form=="feature"){
      cond.en.func <- conditional.config.energy2
    } else {
      stop("Conditional energy formula not properly specified! Use key word feature.function or feature for cond.en.form arguement.")
    }
  }
  for(i in 1:num.states) {

    # Conditional node energy E(Xi | X/Xi)
    node.condtional.energies <- sapply(1:num.nodes,function(xx){
      cond.en.func(
        theta.par                = theta.par,
        config                   = state.space[i,],
        condition.element.number = xx,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)
      # conditional.config.energy(state.space[i,],
      #                           condition.element.number=xx,
      #                           adj.node.list = adjacent.nodes,
      #                           edge.mat = edges.mat,
      #                           one.lgp = node.energies,
      #                           two.lgp = edge.energies,
      #                           ff = ff)
      })

    condtional.energies[i,] <- node.condtional.energies

    #  Complenent conditional node energy E(not-Xi | X/Xi)
    node.complement.condtional.energies <- sapply(1:num.nodes,function(xx){
      cond.en.func(
        theta.par                = theta.par,
        config                   = complement.at.idx(state.space[i,],xx),
        condition.element.number = xx,
        crf                      = crf,
        ff                       = ff,
        printQ                   = FALSE)
      # conditional.config.energy(complement.at.idx(state.space[i,],xx),
      #                           condition.element.number=xx,
      #                           adj.node.list = adjacent.nodes,
      #                           edge.mat = edges.mat,
      #                           one.lgp = node.energies,
      #                           two.lgp = edge.energies,
      #                           ff = ff)
      })

    complement.condtional.energies[i,] <- node.complement.condtional.energies

    # Zs for each conditional:
    node.conditional.Zs <- exp(node.condtional.energies) + exp(node.complement.condtional.energies)
    conditional.Zs[i,]  <- node.conditional.Zs

    # Pr(Xi | X/Xi)
    Prs.ce              <- exp(node.condtional.energies)/node.conditional.Zs
    conditional.Prs[i,] <- Prs.ce

    # Pr(not-Xi | X/Xi)
    #Prs.cce <- exp(node.complement.condtional.energies)/(exp(node.condtional.energies) + exp(node.complement.condtional.energies))
    Prs.cce                        <- 1 - Prs.ce
    complement.conditional.Prs[i,] <- Prs.cce

    # pseudo-liklihood = Prod Pr(Xi | X/Xi)
    pseudo.lik            <- prod(Prs.ce)
    pseudo.likelihoods[i] <- pseudo.lik

  }

  dist.info <- list(
    condtional.energies,
    complement.condtional.energies,
    conditional.Zs,
    conditional.Prs,
    complement.conditional.Prs,
    pseudo.likelihoods
  )
  names(dist.info) <- c(
    "condtional.energies",
    "complement.condtional.energies",
    "conditional.Zs",
    "conditional.Prs",
    "complement.conditional.Prs",
    "pseudo.likelihoods"
  )

  return(dist.info)

}


#' Compute the Kullback-Leibler distance between two distributions
#'
#' Credit: Shamelessly taken from the excellent LaplacesDemon package
#'
#' The function will XXXX
#'
#' @param px vector of probability densities p(x)
#' @param py vector of probability densities p(y)
#'
#' @return The function will XX
#'
#'
#' @export
KLD <- function (px, py, base = exp(1))
{
  if (!is.vector(px))
    px <- as.vector(px)
  if (!is.vector(py))
    py <- as.vector(py)
  n1 <- length(px)
  n2 <- length(py)
  if (!identical(n1, n2))
    stop("px and py must have the same length.")
  if (any(!is.finite(px)) || any(!is.finite(py)))
    stop("px and py must have finite values.")
  if (any(px <= 0))
    px <- exp(px)
  if (any(py <= 0))
    py <- exp(py)
  px[which(px < .Machine$double.xmin)] <- .Machine$double.xmin
  py[which(py < .Machine$double.xmin)] <- .Machine$double.xmin
  px <- px/sum(px)
  py <- py/sum(py)
  KLD.px.py <- px * (log(px, base = base) - log(py, base = base))
  KLD.py.px <- py * (log(py, base = base) - log(px, base = base))
  sum.KLD.px.py <- sum(KLD.px.py)
  sum.KLD.py.px <- sum(KLD.py.px)
  mean.KLD <- (KLD.px.py + KLD.py.px)/2
  mean.sum.KLD <- (sum.KLD.px.py + sum.KLD.py.px)/2
  out <- list(KLD.px.py = KLD.px.py, KLD.py.px = KLD.py.px,
              mean.KLD = mean.KLD, sum.KLD.px.py = sum.KLD.px.py, sum.KLD.py.px = sum.KLD.py.px,
              mean.sum.KLD = mean.sum.KLD, intrinsic.discrepancy = min(sum.KLD.px.py,
                                                                       sum.KLD.py.px))
  return(out)
}
