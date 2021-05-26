#' Port of Schmidt UGM_MRF_makePotentials.m
#'
#' Pretty much a direct port of Schmidt's UGM_MRF_makePotentials.m
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
UGM_MRF_makePotentials<- function(w,crf) {
  # Arg list was originally: (w,nodeMap,edgeMap,edgeStruct)

  nNodes   <- crf$n.nodes
  nEdges   <- crf$n.edges
  edgeEnds <- crf$edges
  nStates  <- crf$n.states
  maxState <- crf$max.state

  nodeMap <- crf$node.par
  edgeMap <- crf$edge.par

  nodePot <- array(0, c(nNodes, maxState))
  for(n in 1:nNodes){
    for(s in 1:nStates[n]){
      if(nodeMap[n,s,1] == 0){
        nodePot[n,s] <- 1
      } else {
        #nodePot(n,s) = exp(w(nodeMap(n,s)));
        nodePot[n,s] <- exp(w[nodeMap[n,s,1]])
      }
    }
  }

  #edgePot = zeros(maxState,maxState,nEdges);
  edgePot <- rep(list(array(0,c(maxState,maxState))),nEdges)
  for(e in 1:nEdges){
    n1 <- edgeEnds[e,1]
    n2 <- edgeEnds[e,1]
    for(es1 in 1:nStates[n1]){
      for(es2 in 1:nStates[n2]){
        if(edgeMap[[e]][es1,es2,1] == 0){
          edgePot[[e]][es1,es2] <- 1
        } else {
          edgePot[[e]][es1,es2] <- exp(w[edgeMap[[e]][es1,es2,1]])
        }
      }
    }
  }

  potentials <- list(nodePot,edgePot)
  names(potentials) <- c("nodePot", "edgePot")

  return(potentials)

}
