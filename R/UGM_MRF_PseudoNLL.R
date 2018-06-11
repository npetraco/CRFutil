#' Experiment: Port of Schmidt UGM_CRF_PseudoNLL.m
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
UGM_MRF_PseudoNLL <- function(w,nInstances,crf) {
  # UGM functionwas: UGM_loss(w,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct)

  nodeMap <- crf$node.par
  edgeMap <- crf$edge.par

  nNodes   <- crf$n.nodes
  #nNodeFeatures = size(Xnode,2);
  #nEdgeFeatures = size(Xedge,2);
  #nEdges = edgeStruct.nEdges;
  nEdges   <- crf$n.edges
  #edgeEnds = edgeStruct.edgeEnds;
  edgeEnds <- crf$edges
  #V = edgeStruct.V;
  #E = edgeStruct.E;
  #nStates = edgeStruct.nStates;
  nStates <- crf$n.states

  NLL <- 0
  g <- numeric(length(w))

  for(i in 1:nInstances){
    # Make potentials
    #nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    mpots <- UGM_MRF_makePotentials(w,crf)
    crf$nodePot <- mpots$nodePot # *************
    crf$edgePot <- mpots$edgePot # *************

    for(n in 1:nNodes){
      # Get Neighbors of node n
      edges <- crf$adj.nodes[[n]]
      print(edges)

      # Compute Probability of Each State with Neighbors Fixed

      for(e in edges){

      }

    }


  }
}
