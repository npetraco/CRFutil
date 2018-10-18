#' XXXX
#'
#' @param XXXX XXXX
#'
#' @return The function will XX
#'
#' @export
symbolic.conditional.energy <- function(config, condition.element.number, crf, ff, format="tex", printQ=FALSE){

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
  if(length(zero.idxs) == 0) { # No zero idxs
    theta.idxs <- param.num.vec
  } else {                     # Drop zero idxs if found
    theta.idxs <- param.num.vec[-which(param.num.vec==0)]
  }
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
    for(i in 2:length(theta.idxs)){
      sne <- paste0(sne," + th_",theta.idxs[i])
      sne.tex <- paste0(sne.tex," + \theta_{",theta.idxs[i],"}")
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
  } else {
    out.eq <- symb.ener
  }

  return(out.eq)

}
