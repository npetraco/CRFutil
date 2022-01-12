#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.exact.nll <- function(par, crf, instances, infer.method = infer.exact, ...){
  mrf.nll(par = par, crf = crf, instances = instances, infer.method = infer.exact, ...)
}


#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.lbp.nll <- function(par, crf, instances, infer.method = infer.lbp, ...){
  mrf.nll(par = par, crf = crf, instances = instances, infer.method = infer.lbp, ...)
}


#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' Work around for CRF bug evaluating nll for different infer.method's other than on chains
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.junction.nll <- function(par, crf, instances, infer.method = infer.junction, ...){
  mrf.nll(par = par, crf = crf, instances = instances, infer.method = infer.junction, ...)
}


#' Decorate initalized mrf-object to make potentials compatible with gRbase
#'
#' Decorate initalized mrf-object to make potentials compatible with gRbase
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.gRbase.potentials <- function(crf, node.names, state.nmes=NULL){

  if(is.null(state.nmes)) {
    state.nmes.loc <- paste0("s",1:crf$max.state)
  } else {
    state.nmes.loc <- state.nmes
  }

  # Decorate node potentials:
  gRbase.node.potentials     <- rep(list(NULL),crf$n.nodes) #node Psi's
  gRbase.node.log.potentials <- rep(list(NULL),crf$n.nodes) #node psi's (one-body energies)
  for(i in 1:crf$n.nodes){
    node.levs                       <- list(state.nmes.loc)
    names(node.levs)                <- node.names[i]
    # train.mrf puts an extra dimension onto node.pot. To account for that:
    if(length(dim(crf$node.pot))==3) {
      vls <- crf$node.pot[i,,1]
    } else {
      vls <- crf$node.pot[i,]
    }
    gRbase.node.potentials[[i]]     <- ar_new(node.names[i], levels=node.levs, values=c(vls))
    gRbase.node.log.potentials[[i]] <- log(gRbase.node.potentials[[i]])
  }

  # Decorate edge potentials:
  gRbase.edge.potentials     <- rep(list(NULL),crf$n.edges) # edge Psi's
  gRbase.edge.log.potentials <- rep(list(NULL),crf$n.edges) # edge psi's (two-body energies)
  for(i in 1:crf$n.edges){
    e1                              <- node.names[crf$edges[i,1]]
    e2                              <- node.names[crf$edges[i,2]]
    node.levs                       <- list(state.nmes.loc,state.nmes.loc)
    names(node.levs)                <- c(e1,e2)
    gRbase.edge.potentials[[i]]     <- ar_new(c(e1,e2), levels=node.levs, values=as.numeric(crf$edge.pot[[i]]))
    gRbase.edge.log.potentials[[i]] <- log(gRbase.edge.potentials[[i]])
  }

  potential.info        <- list(gRbase.node.potentials,
                                gRbase.edge.potentials,
                                gRbase.node.log.potentials,
                                gRbase.edge.log.potentials)
  names(potential.info) <- c("node.potentials","edge.potentials","node.energies","edge.energies")

  return(potential.info)

}


#' Decorate fit node and edge marginal beliefs to make compatible with gRbase
#'
#' Decorate fit node and edge marginal beliefs to make compatible with gRbase
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.gRbase.beliefs <- function(inference.obj, node.names, edge.mat, state.nmes=NULL){

  num.nodes  <- nrow(inference.obj$node.bel)
  num.edges  <- length(inference.obj$edge.bel)
  num.states <- ncol(inference.obj$node.bel)

  if(is.null(state.nmes)) {
    state.nmes.loc <- paste0("s",1:num.states)
  } else {
    state.nmes.loc <- state.nmes
  }

  # Decorate node beliefs:
  gRbase.node.bels <- rep(list(NULL), num.nodes) #node beliefs
  for(i in 1:num.nodes){
    node.levs                 <- list(state.nmes.loc)
    names(node.levs)          <- node.names[i]
    gRbase.node.bels[[i]]     <- ar_new(node.names[i], levels=node.levs, values=c(inference.obj$node.bel[i,]))
  }

  # Decorate edge beliefs:
  gRbase.edge.bels <- rep(list(NULL), num.edges) # edge beliefs
  for(i in 1:num.edges){
    e1                        <- node.names[edge.mat[i,1]]
    e2                        <- node.names[edge.mat[i,2]]
    node.levs                 <- list(state.nmes.loc,state.nmes.loc)
    names(node.levs)          <- c(e1,e2)
    gRbase.edge.bels[[i]]     <- ar_new(c(e1,e2), levels=node.levs, values=as.numeric(inference.obj$edge.bel[[i]]))
  }

  belief.info        <- list(gRbase.node.bels,
                             gRbase.edge.bels)
  names(belief.info) <- c("node.beliefs","edge.beliefs")

  return(belief.info)

}

#' Instantiate an empty field
#'
#' XXXX
#'
#' XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.empty.field <- function(graph.eq=NULL, adj.mat=NULL, parameterization.typ="standard", node.par=NULL, edge.par=NULL, plotQ=FALSE) {

  if(is.null(graph.eq) & is.null(adj.mat)){
    stop("Specify a model!")
  }

  if(!is.null(graph.eq) & !is.null(adj.mat)){
    stop("Specify model as either a graph eq OR an adjacency matrix.")
  }

  # if(num.states != 2) {
  #   warning("Caution! Number of node states not equal to 2. Most functionality in CRFutil assumes there are only 2 states!")
  # }
  num.states <- 2 # Just assume two states per node for now. Parameterization code would have to be much more elaborate otherwise.

  if(!is.null(graph.eq)){
    adjm <- ug(graph.eq, result="matrix")
    # Check the graph:
    #gph <- ug(grphf, result = "graph")

  } else {
    adjm <- adj.mat
  }
  new.crf <- make.crf(adjm, num.states)
  new.crf <- make.features(new.crf)

  if(parameterization.typ == "standard") {

    # One parameter per node, one parameter per edge
    num.pars <- new.crf$n.nodes + new.crf$n.edges
    new.crf <- make.par(new.crf, num.pars)

    par.count <- 1
    for(i in 1:new.crf$n.nodes) {
      new.crf$node.par[par.count,1,] <- par.count
      par.count <- par.count + 1
    }
    for(i in 1:new.crf$n.edges){
      new.crf$edge.par[[i]][1,1,1] <- par.count
      new.crf$edge.par[[i]][2,2,1] <- par.count
      par.count <- par.count + 1
    }

  } else if(parameterization.typ == "flexible") {

    # One parameter per node, two parameters per edge
    num.pars <- new.crf$n.nodes + 2*new.crf$n.edges
    new.crf <- make.par(new.crf, num.pars)

    par.count <- 1
    for(i in 1:new.crf$n.nodes) {
      new.crf$node.par[par.count,1,] <- par.count
      par.count <- par.count + 1
    }
    for(i in 1:new.crf$n.edges){
      new.crf$edge.par[[i]][1,1,1] <- par.count
      par.count <- par.count + 1
      new.crf$edge.par[[i]][2,2,1] <- par.count
      par.count <- par.count + 1
    }


  } else if(parameterization.typ == "ising1") {

    # No node parameters and one parameter for all the edges
    num.pars <- 1
    new.crf <- make.par(new.crf, num.pars)

    for(i in 1:new.crf$n.edges){
      new.crf$edge.par[[i]][1,1,1] <- 1
      new.crf$edge.par[[i]][2,2,1] <- 1
    }


  } else if(parameterization.typ == "ising2") {

    # One parameter for all the nodes, one parameter for all the edges
    num.pars <- 2
    new.crf <- make.par(new.crf, num.pars)

    for(i in 1:new.crf$n.nodes) {
      new.crf$node.par[i,1,] <- 1
    }
    for(i in 1:new.crf$n.edges){
      new.crf$edge.par[[i]][1,1,1] <- 2
      new.crf$edge.par[[i]][2,2,1] <- 2
    }

  } else if(parameterization.typ == "ising3") {

    # Different parameters for all the nodes, one parameter for all the edges
    num.pars <- new.crf$n.nodes+1
    new.crf  <- make.par(new.crf, num.pars)

    for(i in 1:new.crf$n.nodes) {
      new.crf$node.par[i,1,] <- i
    }
    for(i in 1:new.crf$n.edges){
      new.crf$edge.par[[i]][1,1,1] <- new.crf$n.nodes+1
      new.crf$edge.par[[i]][2,2,1] <- new.crf$n.nodes+1
    }


  } else if(parameterization.typ == "custom") {

    if(is.null(node.par) & is.null(edge.par)){
      stop("Custom parameterization specified but no node or edge pars given!")
    }

    num.pars <- max(c(as.numeric(node.par), unlist(edge.par)))
    new.crf <- make.par(new.crf, num.pars)

    if(!is.null(node.par)){
      new.crf$node.par <- node.par
    }
    if(!is.null(edge.par)){
      new.crf$edge.par <- edge.par
    }

  } else {
    stop("Specify parameterization choice: standard, flexible, ising1, ising2, ising3, or custom!")
  }

  #dump.crf(new.crf)
  if(plotQ==TRUE){
    new.crf.gp <- as(adjm,"graphNEL")
    if(!is.null(dev.list())){
      dev.off()
    }
    iplot(new.crf.gp)
  }

  return(new.crf)

}

#' (Deep) Copy and return a new independent crf object
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
copy.crf <- function(crf, plotQ=FALSE){

  adj.mat <- array(0, c(crf$n.nodes, crf$n.nodes))
  for(i in 1:length(crf$adj.nodes)) {
    idx <- crf$adj.nodes[[i]]
    adj.mat[idx[1], idx[2]] <- 1
  }
  adj.mat <- t(adj.mat) + adj.mat # symmetrize assuming start is only upper/lower triangle
  colnames(adj.mat) <- 1:crf$n.nodes
  rownames(adj.mat) <- 1:crf$n.nodes

  if(max(adj.mat) > 1) {
    print(adj.mat)
    stop("Something went wrong with reconstructing the adjacency matrix!")
  }

  new.crf <- make.crf(adj.mat, crf$n.states)
  new.crf <- make.features(new.crf)
  new.crf <- make.par(new.crf, crf$n.par)
  crf.attrib.nms <- names(crf)
  for(i in 1:length(crf.attrib.nms)){
    new.crf[[crf.attrib.nms[i]]] <- crf[[crf.attrib.nms[i]]]
  }

  if(plotQ==TRUE){
    new.crf.gp <- as(adj.mat,"graphNEL")
    if(!is.null(dev.list())){
      dev.off()
    }
    iplot(new.crf.gp)
  }

  return(new.crf)

}


#' Print out all contents of a crf object
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
dump.crf <- function(crf){

  crf.attrib.nms <- names(crf)

  for(i in 1:length(crf.attrib.nms)){
    print("----------------")
    print(crf.attrib.nms[i])
    print(crf[[crf.attrib.nms[i]]])
  }

}
