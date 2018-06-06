#' Experiment Port of Schmidt UGM_MRF_computeSuffStat.m
#'
#' Experiment
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
MRF_Stat_expt <- function(crf, instances) {

  nNodes   <- crf$n.nodes
  nEdges   <- crf$n.edges
  edgeEnds <- crf$edges
  nStates  <- crf$n.states
  nParams  <- crf$n.par

  nInstances <- nrow(instances)
  suffStat   <- numeric(nParams)

  nodeMap <- crf$node.par
  edgeMap <- crf$edge.par
  for(i in 1:nInstances){
    y <- instances[i,]
    for(n in 1:nNodes){
      if(nodeMap[n,y[n],1] > 0){
        suffStat[nodeMap[n,y[n],1]] <- suffStat[nodeMap[n,y[n],1]] + 1
      }
    }
    for(e in 1:nEdges){
      n1 <- edgeEnds[e,1]
      n2 <- edgeEnds[e,2]
      if(edgeMap[[e]][y[n1],y[n2],1] > 0) {
        suffStat[edgeMap[[e]][y[n1],y[n2],1]] <- suffStat[edgeMap[[e]][y[n1],y[n2],1]] + 1
      }
    }
  }

  return(suffStat)
}
