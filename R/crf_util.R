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
    state.nmes <- paste0("s",1:crf$max.state)
  }

  # Decorate node potentials:
  gRbase.node.potentials     <- rep(list(NULL),crf$n.nodes) #node Psi's
  gRbase.node.log.potentials <- rep(list(NULL),crf$n.nodes) #node psi's (one-body energies)
  for(i in 1:crf$n.nodes){
    node.levs                       <- list(state.nmes)
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
    node.levs                       <- list(state.nmes,state.nmes)
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
    state.nmes <- paste0("s",1:num.states)
  }

  # Decorate node beliefs:
  gRbase.node.bels <- rep(list(NULL), num.nodes) #node beliefs
  for(i in 1:num.nodes){
    node.levs                 <- list(state.nmes)
    names(node.levs)          <- node.names[i]
    gRbase.node.bels[[i]]     <- ar_new(node.names[i], levels=node.levs, values=c(inference.obj$node.bel[i,]))
  }

  # Decorate edge beliefs:
  gRbase.edge.bels <- rep(list(NULL), num.edges) # edge beliefs
  for(i in 1:num.edges){
    e1                        <- node.names[edge.mat[i,1]]
    e2                        <- node.names[edge.mat[i,2]]
    node.levs                 <- list(state.nmes,state.nmes)
    names(node.levs)          <- c(e1,e2)
    gRbase.edge.bels[[i]]     <- ar_new(c(e1,e2), levels=node.levs, values=as.numeric(inference.obj$edge.bel[[i]]))
  }

  belief.info        <- list(gRbase.node.bels,
                             gRbase.edge.bels)
  names(belief.info) <- c("node.beliefs","edge.beliefs")

  return(belief.info)

}

