#' Experiment: Port of Schmidt UGM_MRF_computeSuffStat.m
#'
#' Experiment: The sufficient statistics seem to be the number
#' of non-zero energy parameters... See loops in function below
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
UGM_MRF_computeSuffStat <- function(crf, instances) {
  # Arg list was originally: (Y,nodeMap,edgeMap,edgeStruct)
  # Changed to match CRF function mrf.stat arg list

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


#' Experiment: Port of Schmidt UGM_MRF_makePotentials.m
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


#' Experiment: Port of Schmidt UGM_MRF_NLL.m
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
UGM_MRF_NLL <- function(w,nInstances,suffStat,crf,inferFunc) {
  # CRF function: mrf.nll(par, crf, instances, infer.method = infer.chain, ...)
  # UGM function arg list was: w,nInstances,suffStat,nodeMap,edgeMap,edgeStruct,inferFunc

  nNodes   <- crf$n.nodes
  nEdges   <- crf$n.edges
  edgeEnds <- crf$edges
  nStates  <- crf$n.states

  #w       <- crf$par
  nodeMap <- crf$node.par
  edgeMap <- crf$edge.par

  # Make potentials
  #[nodePot,edgePot] = UGM_MRF_makePotentials(w,nodeMap,edgeMap,edgeStruct);
  mpots <- UGM_MRF_makePotentials(w, crf)
  #nodePot <- mpots$nodePot
  #edgePot <- mpots$edgePot
  crf$nodePot <- mpots$nodePot
  crf$edgePot <- mpots$edgePot

  # Compute marginals and logZ
  #[nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
  infer.info <- inferFunc(crf)
  logZZ <- infer.info$logZ
  print(infer.info)

  # Compute NLL
  #NLL = -w'*suffStat + nInstances*logZ;
  NLL <- -w %*% suffStat + nInstances*logZZ
  print(NLL)


}
