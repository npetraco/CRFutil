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
  crf$nodePot <- mpots$nodePot # Has to be done here because of arg structure to inferFunc
  crf$edgePot <- mpots$edgePot # Has to be done here because of arg structure to inferFunc

  # Compute marginals and logZ
  #[nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
  infer.info <- inferFunc(crf)
  logZZ   <- infer.info$logZ
  nodeBel <- infer.info$node.bel
  edgeBel <- infer.info$edge.bel

  # Compute negative log liklihood:
  #NLL = -w'*suffStat + nInstances*logZ;
  NLL <- -w %*% suffStat + nInstances*logZZ

  # Gradient of NLL:
  g <- -suffStat
  for(n in 1:nNodes){
    for(s in 1:nStates[n]){
      if(nodeMap[n,s,1] > 0){
        #g(nodeMap(n,s)) = g(nodeMap(n,s)) + nInstances*nodeBel(n,s);
        g[nodeMap[n,s,1]] <- g[nodeMap[n,s,1]] + nInstances*nodeBel[n,s]
      }
    }
  }
  for(e in 1:nEdges){
    n1 <- edgeEnds[e,1]
    n2 <- edgeEnds[e,2]
    for(es1 in 1:nStates[n1]){
      for(es2 in 1:nStates[n2]){
        if(edgeMap[[e]][es1,es2,1] > 0){
          #g(edgeMap(s1,s2,e)) = g(edgeMap(s1,s2,e)) + nInstances*edgeBel(s1,s2,e);
          g[edgeMap[[e]][es1,es2,1]] <- g[edgeMap[[e]][es1,es2,1]] + nInstances*edgeBel[[e]][es1,es2]
        }
      }
    }
  }

  nll.info <- list(NLL,g)
  names(nll.info) <- c("neg.log.lik","grad")

  return(nll.info)

}
