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
#' Compute total energy of a configuration. Assumes node/edge energies are in gRbase format   *****GENERALIZE
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

#' Energy function for energy of a conditional configuration of states E(Xi | X/Xi).
#' Uses the feature function formulation. Eliminates dependency on gRbase format for potentials
#'
#' Assumes log-potentials are in gRbase format
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
conditional.config.energy <- function(theta.par=NULL, config, condition.element.number, crf, ff, printQ=FALSE) {

  # Use input theta if one is supplied, otherwise use theta that is in crf object
  if(is.null(theta.par)) {
    node.pot      <- crf$node.pot
    edge.pot      <- crf$edge.pot
  } else {
    pots <- make.pots(parms=theta.par, crf = crf, rescaleQ=F, replaceQ=F, format = "regular", printQ=F)
    node.pot <- pots[[1]]
    edge.pot <- pots[[2]]
  }

  adj.node.list <- crf$adj.nodes
  edge.mat      <- crf$edges

  # One-body energy (log node-potentials)
  e.one <- Eone(config[condition.element.number], log(node.pot[condition.element.number,]), ff)

  # Nodes connected to the condition node Xi
  adj.nodes <- adj.node.list[[condition.element.number]]

  e.two <- 0
  for(i in 1:length(adj.nodes)){
    edg         <- sort(c(condition.element.number, adj.nodes[i]))
    edg.pot.idx <- row.match(edg, table = edge.mat)
    e.two       <- e.two + Etwo(config[edg[1]], config[edg[2]], log(edge.pot[[edg.pot.idx]]), ff)

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

#' Older implementation of Energy function for energy of a conditional configuration of states E(Xi | X/Xi).
#' Uses the feature function formulation.
#'
#' Assumes log-potentials are in gRbase format
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
conditional.config.energy.old <- function(config, condition.element.number, adj.node.list, edge.mat, one.lgp, two.lgp, ff, printQ=FALSE) {

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


#' Energy function for energy of a conditional configuration of states E(Xi | X/Xi).
#' Uses the feature formulation.
#'
#' Assumes log-potentials are in gRbase format
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
conditional.config.energy2 <- function(theta.par=NULL, config, condition.element.number, crf, ff, printQ=FALSE) {

  # Node component:
  l     <- get.par.idx(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)
  phi.l <- phi.component(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)

  #print(l)
  #print(l.phi)
  if(l==0){
    if(phi.l != 0) {
      stop("Something is wrong. l-idx is 0 but its phi is not.")
    }
  }

  if(phi.l==0){
    if(l != 0) {
      stop("Something is wrong. phi for this l is 0 but l is not.")
    }
  }

  # Use input theta if one is supplied, otherwise use theta that is in crf object
  if(is.null(theta.par)) {
    local.par <- crf$par
  } else {
    local.par <- theta.par
  }

  configE <- 0
  if(phi.l !=0) {
    configE <- configE + local.par[l]
  }

  adj.nodes <- crf$adj.nodes[[condition.element.number]]
  for(ii in 1:length(adj.nodes)) {

    edge.nods <- sort(c(condition.element.number, adj.nodes[ii]))

    if(printQ==T){
      # print(config)
      print(paste("ii: ", ii))
      print("Edge-Nodes:")
      print(edge.nods[1])
      print(edge.nods[2])
      # print(crf$edge.par)
      # print(crf$edges)
    }

    k <- get.par.idx(config = config,
                     i        = edge.nods[1],
                     j        = edge.nods[2],
                     edge.par = crf$edge.par,
                     edge.mat = crf$edges,
                     ff       = ff)
    if(printQ==T){
      print("phi.component calls get.par.idx again:")
    }
    phi.k <- phi.component(config = config,
                     i        = edge.nods[1],
                     j        = edge.nods[2],
                     edge.par = crf$edge.par,
                     edge.mat = crf$edges,
                     ff       = ff)

    if(k==0){
      if(phi.k != 0) {
        stop(paste("Something is wrong. k-idx is 0 but its phi is not. k=", k, " phi.k=", phi.k))
      }
    }

    if(phi.k==0){
      if(k != 0) {
        stop(paste("Something is wrong. phi for this k is 0 but k is not. k=", k, " phi.k=", phi.k))
      }
    }

    # print(k)
    # print(phi.k)

    if(phi.k !=0) {
      configE <- configE + local.par[k]
    }

  }

 return(configE)

}


#' Spit out latex expression for E(Xi|X/Xi) given i, X and a crf. Basically used feature formulation for conditional energy.
#'
#' Can also spit out a condtional phi vector  (i.e. and alpha vector) if format="conditional.phi" is specified
#'
#' @param XXXX XXXX
#'
#' @return The function will XX
#'
#' @export
symbolic.conditional.energy <- function(config, condition.element.number, crf, ff, format="tex", printQ=FALSE){    # **** NEEDS TO BE C

  param.num.vec <- NULL
  cond.phi.vec  <- NULL

  # Parameter (if any) associated with conditioned node
  l     <- get.par.idx(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)
  phi.l <- phi.component(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)
  if(printQ==TRUE){
    print(paste0("For node i: ", condition.element.number, " in state Xi=", config[condition.element.number], ", param# assoc l=", l, " and thus phi_l=", phi.l))
  }

  param.num.vec <- c(param.num.vec, l)
  cond.phi.vec  <- c(cond.phi.vec, phi.l)

  # Parameter (if any) associated with edge containing conditioned node
  adj.nodes <- crf$adj.nodes[[condition.element.number]]
  if(printQ==TRUE){
    print(paste0("The nodes below are connected to node i=", condition.element.number))
    print(adj.nodes)
  }

  for(ii in 1:length(adj.nodes)) {

    edge.nods <- sort(c(condition.element.number, adj.nodes[ii]))

    k <- get.par.idx(
      config   = config,
      i        = edge.nods[1],
      j        = edge.nods[2],
      edge.par = crf$edge.par,
      edge.mat = crf$edges,
      ff       = ff)

    phi.k <- phi.component(
      config   = config,
      i        = edge.nods[1],
      j        = edge.nods[2],
      edge.par = crf$edge.par,
      edge.mat = crf$edges,
      ff       = ff)

    param.num.vec <- c(param.num.vec, k)
    cond.phi.vec  <- c(cond.phi.vec, phi.k)

    if(printQ==TRUE){
      print(paste0("For edge #", ii, " edge: ", edge.nods[1],"-",edge.nods[2],
                   " in states Xi=", config[edge.nods[1]], " Xj=", config[edge.nods[2]],
                   ", param# assoc k=", k, " and thus phi_k=", phi.k))
    }
  }

  if(printQ==TRUE){
    print("Parameter number vec: ")
    print(param.num.vec)
  }

  zero.idxs <- which(param.num.vec==0)
  #print(zero.idxs)
  if(length(zero.idxs) == 0) { # No zero idxs
    theta.idxs <- param.num.vec
  } else {                     # Drop zero idxs if found
    theta.idxs <- param.num.vec[-which(param.num.vec==0)]
  }
  #print(theta.idxs)

  if(printQ==TRUE){
    print("Final theta indices: ")
    print(theta.idxs)
  }

  if(length(theta.idxs) == 0) {
    sne <- "0"
    sne.tex <- "0"
  } else {
    sne <- paste0("th_",theta.idxs[1])
    sne.tex <- paste0("\theta_{",theta.idxs[1],"}")
    if(length(theta.idxs) >=2 ) {
      for(i in 2:length(theta.idxs)){
        sne <- paste0(sne," + th_",theta.idxs[i])
        sne.tex <- paste0(sne.tex," + \theta_{",theta.idxs[i],"}")
      }
    }
  }

  #sne <- paste0("th_",theta.idxs)
  #print(sne)
  # do all drop case
  symb.ener     <- paste0("E(X_",condition.element.number,"|X/X_",condition.element.number,") = ", sne)
  symb.ener.tex <- paste0("E(X_{",condition.element.number,"}|{\bf X}/X_{",condition.element.number,"}) = ", sne.tex)

  if(printQ==TRUE){
    print(symb.ener)
    print(symb.ener.tex)
  }

  if(format=="tex"){
    out.eq <- symb.ener.tex
  } else if(format=="conditional.phi"){
    conditional.phi.tmp <- array(0,c(1,crf$n.par))
    conditional.phi.tmp[theta.idxs] <- 1
    out.eq <- conditional.phi.tmp
    colnames(out.eq) <- 1:crf$n.par
    #print(theta.idxs)
    #print(cond.phi.vec)

  } else {
    out.eq <- symb.ener
  }

  return(out.eq)

}


#' Gradient of a condtional energy for a node and a configuration
#'
#' @param XXXX XXXX
#' @details This function compute the vector: \nabla_{\boldsymbol \theta}E(X_i|{\bf X}\slash X_i). The vector is the same length as the
#'          number of parameters.
#' @return The function will XX
#'
#' @export
conditional.energy.gradient <- function(config, condition.element.number, crf, ff, printQ=FALSE){

  param.num.vec <- NULL
  cond.phi.vec  <- NULL

  # Parameter (if any) associated with conditioned node
  l     <- get.par.idx(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)
  phi.l <- phi.component(config = config, i=condition.element.number, node.par=crf$node.par, ff=ff)
  if(printQ==TRUE){
    print(paste0("For node i: ", condition.element.number, " in state Xi=", config[condition.element.number], ", param# assoc l=", l, " and thus phi_l=", phi.l))
  }

  param.num.vec <- c(param.num.vec, l)
  cond.phi.vec  <- c(cond.phi.vec, phi.l)

  # Parameter (if any) associated with edge containing conditioned node
  adj.nodes <- crf$adj.nodes[[condition.element.number]]
  if(printQ==TRUE){
    print(paste0("The nodes below are connected to node i=", condition.element.number))
    print(adj.nodes)
  }

  for(ii in 1:length(adj.nodes)) {

    edge.nods <- sort(c(condition.element.number, adj.nodes[ii]))

    k <- get.par.idx(
      config   = config,
      i        = edge.nods[1],
      j        = edge.nods[2],
      edge.par = crf$edge.par,
      edge.mat = crf$edges,
      ff       = ff)

    phi.k <- phi.component(
      config   = config,
      i        = edge.nods[1],
      j        = edge.nods[2],
      edge.par = crf$edge.par,
      edge.mat = crf$edges,
      ff       = ff)

    param.num.vec <- c(param.num.vec, k)
    cond.phi.vec  <- c(cond.phi.vec, phi.k)

    if(printQ==TRUE){
      print(paste0("For edge #", ii, " edge: ", edge.nods[1],"-",edge.nods[2],
                   " in states Xi=", config[edge.nods[1]], " Xj=", config[edge.nods[2]],
                   ", param# assoc k=", k, " and thus phi_k=", phi.k))
    }
  }

  conditional.grad <- numeric(crf$n.par) # The components of this are \alpha_{k_{[ \sim i]}}({\bf X})
  theta.counts <- table(param.num.vec)
  thetas <- as.numeric(names(theta.counts))
  for(i in 1:length(thetas)){
    conditional.grad[ thetas[i] ] <- theta.counts[i]
  }

  param.info        <- list(cond.phi.vec, param.num.vec, conditional.grad)
  names(param.info) <- c("conditional.phi.components", "parameter.numbers", "conditional.grad")
  return(param.info)

}
