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

  prod.edge.pots <- tableMult(gRbase.edge.potentials[[2]], gRbase.edge.potentials[[1]])
  if(num.edges > 2){
    for(i in 3:num.edges){
      prod.edge.pots <- tableMult(prod.edge.pots, gRbase.edge.potentials[[i]])
    }
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


#' Compute the pseudolikelihoods from the node and edge energies
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
pseudolikelihoods.from.energies <- function(state.space, adjacent.nodes, edges.mat, node.energies, edge.energies, conditional.energy.func, ff){


  num.nodes       <- ncol(state.space)
  num.states      <- nrow(state.space)

  condtional.energies
  complement.condtional.energies
  conditional.Prs
  complement.conditional.Prs
  pseudo.likelihoods
  for(i in 1:num.states) {

    #E(Xi | X/Xi)
    node.condtional.energies <- sapply(1:num.nodes,function(xx){
      conditional.config.energy(state.space[i,],
                                condition.element.number=xx,
                                adj.node.list = adjacent.nodes,
                                edge.mat = edges.mat,
                                one.lgp = node.energies,
                                two.lgp = edge.energies,
                                ff = ff)})

    #E(not-Xi | X/Xi)
    node.complement.condtional.energies <- sapply(1:num.nodes,function(xx){
      conditional.config.energy(complement.at.idx(state.space[i,],xx),
                                condition.element.number=xx,
                                adj.node.list = adjacent.nodes,
                                edge.mat = edges.mat,
                                one.lgp = node.energies,
                                two.lgp = edge.energies,
                                ff = ff)})

    Prs.ce  <- exp(node.condtional.energies)/(exp(node.condtional.energies) + exp(cce))
    Prs.cce <- exp(node.complement.condtional.energies)/(exp(node.condtional.energies) + exp(node.complement.condtional.energies))

    pseudo.lik <- prod(Prs.ce)

  }
  #state.energies  <- sapply(1:num.states, function(xx){energy.func(state.space[xx,], edges.mat, node.energies, edge.energies, ff)})
  #logZZ           <- logsumexp2(state.energies)
  #log.state.probs <- state.energies-logZZ

  #dist.info <- list(exp(log.state.probs), logZZ)
  #names(dist.info) <- c("state.probs", "logZ")
  return(dist.info)

}
