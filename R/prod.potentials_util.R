#' Product of potentials function over a configuration of states
#'
#' Product of potentials function. Expects potentials in gRbase format
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
config.potential <- function(config, edges.mat, node.pots, edge.pots) {

  num.nodes <- length(config)
  num.edges <- nrow(edges.mat)

  # Product of realized configuration node potentials
  prod.node.pots <- 1
  for(i in 1:num.nodes){
    prod.node.pots <- prod.node.pots * node.pots[[i]][config[i]]
  }

  # Product of realized configuration edge potentials
  prod.edge.pots <- 1
  for(i in 1:num.edges){
    prod.edge.pots <- prod.edge.pots * edge.pots[[i]][config[edges.mat[i,1]], config[edges.mat[i,2]]]
  }

  config.prod.pot <- as.numeric(prod.node.pots*prod.edge.pots)

  return(config.prod.pot)
}
