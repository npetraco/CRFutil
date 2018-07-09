#' One-body (node) energy function
#'
#' One-body (node) wnergy function
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
Eone  <- function(yA, tA, ff){ tA %*% ff(yA) }


#' Two-body (edge) energy function
#'
#' Two-body (edge) energy function
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
Etwo  <- function(yA, yB, wAB, ff){ ff(yA) %*% wAB %*% ff(yB) }


#' Energy function for a configuration of states
#'
#' Compute total energy of a configuration. Assumes node/edge energies are in gRbase format
#'
#' The function will XXXX
#'
#' @param config    A node configuration (state) vector
#' @param edges.mat Matrix of connected node edges
#' @param two.lgp   Log node potentials (one-body energies)
#' @param two.lgp   Log edge potentials (two-body energies)
#' @param ff        The feature function
#' @return The function will XX
#'
#'
#' @export
config.energy <- function(config, edges.mat, one.lgp, two.lgp, ff) {

  num.nodes <- length(config)
  num.edges <- nrow(edges.mat)

  # Sum One-body energies (log node-potentials)
  e.one <- 0
  for(i in 1:num.nodes){
    e.one <- e.one + Eone(config[i], one.lgp[[i]], ff)
  }

  # Sum Two-body energies (log edge-potentials)
  e.two <- 0
  for(i in 1:num.edges){
    e.two <- e.two + Etwo(config[edges.mat[i,1]], config[edges.mat[i,2]], two.lgp[[i]], ff)
  }

  ener <- as.numeric(e.one + e.two)

  return(ener)
}

#' Energy function for energy of a conditional configuration of states E(Xi | X/Xi)
#'
#' Compute
#'
#' The function will XXXX
#'
#' @param config    A node configuration (state) vector
#' @param condition.element.number i of E(Xi | X/Xi)
#' @param adj.node.list XXXX
#' @param edges.mat Matrix of connected node edges
#' @param two.lgp   Log node potentials (one-body energies)
#' @param two.lgp   Log edge potentials (two-body energies)
#' @param ff        The feature function
#' @return The function will XX
#'
#'
#' @export
conditional.config.energy <- function(config, condition.element.number, adj.node.list, edge.mat, one.lgp, two.lgp, ff, printQ=FALSE) {

  #num.nodes <- length(config)

  # One-body energy (log node-potentials)
  e.one <- Eone(config[condition.element.number], one.lgp[[condition.element.number]], ff)

  # Nodes connected to the condition node Xi
  adj.nodes <- adj.node.list[[condition.element.number]]

  e.two <- 0
  for(i in 1:length(adj.nodes)){
    edg         <- sort(c(condition.element.number, adj.nodes[i]))
    edg.pot.idx <- row.match(edg, table = edge.mat)
    e.two       <- e.two + Etwo(config[edg[1]], config[edg[2]], two.lgp[[edg.pot.idx]], ff)

    # For checks:
    if(printQ==TRUE){
      print(config)
      print(paste("CE-idx:", condition.element.number))
      print(paste("Edge 1:", edg[1], "    Edge 2:", edg[2]))
      print(paste("St-Edge 1:", config[edg[1]], "   St-Edge 2:", config[edg[2]]))
      print(paste("Edge Pot #:", edg.pot.idx))
      print(paste("e.two:",e.two))
      print(paste("e.one:",e.one))
      print("----------------------")
    }

  }

  ener <- as.numeric(e.one + e.two)
  return(ener)

}


#' Energy function for energy of a conditional configuration of states E(Xi | X/Xi), computed in terms of phi features
#'
#' Compute
#'
#' The function will XXXX
#'
#' @param config    A node configuration (state) vector
#' @param condition.element.number i of E(Xi | X/Xi)
#' @param crf       CRF object
#' @return          The function will XX
#'
#'
#' @export
conditional.config.energy2 <- function(config, phi.config = NULL, condition.element.number, crf, ff = NULL) {

  if(is.null(crf$nodes2pars)){
    stop("Compute nodes2parameters list and store in crf object!")
  }

  # Compute phi features of config if they were not passed in:
  if(is.null(phi.config)){

    if(is.null(ff)){
      stop("Include a feature function (ff) or a phi.config!")
    }
    phi.config <- phi.features(
      config    = config,
      edges.mat = crf$edges,
      node.par  = crf$node.par,
      edge.par  = crf$edge.par,
      ff        = ff
    )
  }

  ener <- crf$par[ crf$nodes2pars[[condition.element.number]] ] %*% phi.config[ crf$nodes2pars[[condition.element.number]] ]

  return(ener)

}


#' Energy function for energy of a conditional configuration of states E(Xi | X/Xi), computed in terms of phi features
#' This time however pass in a par (theta) vector instead of using the one with the crf object/
#'
#' The function will XXXX
#'
#' @param config    A node configuration (state) vector
#' @param condition.element.number i of E(Xi | X/Xi)
#' @param crf       CRF object
#' @return          The function will XX
#'
#'
#' @export
conditional.config.energy2a <- function(par, config, phi.config = NULL, condition.element.number, crf, ff = NULL) {

  if(is.null(crf$nodes2pars)){
    stop("Compute nodes2parameters list and store in crf object!")
  }

  # Compute phi features of config if they were not passed in:
  if(is.null(phi.config)){

    if(is.null(ff)){
      stop("Include a feature function (ff) or a phi.config!")
    }
    phi.config <- phi.features(
      config    = config,
      edges.mat = crf$edges,
      node.par  = crf$node.par,
      edge.par  = crf$edge.par,
      ff        = ff
    )
  }

  ener <- par[ crf$nodes2pars[[condition.element.number]] ] %*% phi.config[ crf$nodes2pars[[condition.element.number]] ]

  return(ener)

}
