#' Convert loglin, loglm or dModel (all really loglin) object to a crf object
#'
#' Convert loglin, loglm or dModel (all really loglin) object to a crf object
#'
#' @details Do glm, glmnet /logistic/poisson model conversion in a different function
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
loglin2crf <- function(log.linear.model.obj, standard.potentialsQ=T, plotQ=F) {

  model.class <- class(log.linear.model.obj)[[1]]

  if(model.class == "list") {          # For loglin() model
    model.coefs <- log.linear.model.obj$param
  } else if(model.class  == "loglm") { # For loglm() model
    model.coefs <- log.linear.model.obj$param
  } else if(model.class == "dModel") { # For gRim dmod() model
    model.coefs <- log.linear.model.obj$fitinfo$param
  } else {
    stop("Model type must be list (loglin), loglm (MASS) or dModel (gRim)!")
  }

  # Parse-out which elements in the model.coefs list are for nodes and which are for edges:
  node.and.edge.names <- names(model.coefs)

  # Indicator vector of whether parameters are for a node or an edge.
  # FALSE means node, TRUE means edge
  edgeQ <- sapply(1:length(node.and.edge.names), function(xx){grepl( ".", node.and.edge.names[xx], fixed = TRUE)})

  edge.idxs     <- which(edgeQ == TRUE)
  intercept.idx <- which(node.and.edge.names == "(Intercept)") # May or may not be an intercept
  node.idxs     <- 1:length(node.and.edge.names)
  node.idxs     <- node.idxs[-c(intercept.idx, edge.idxs)]

  # Build edge matrix from input log linear object.
  # **NOTE: Number by canonical node order extracted from names(model.coefs)
  node.names     <- node.and.edge.names[node.idxs]
  node.names.mat <- data.frame(1:length(node.names), node.names)
  colnames(node.names.mat) <- c("node.idx", "node.name")

  edge.mat.loc       <- array(NaN, c(length(edge.idxs), 2))
  edge.mat.names.loc <- array("", c(length(edge.idxs), 2))
  for(i in 1:length(edge.idxs)){
    edg.loc <- node.and.edge.names[edge.idxs[i]]         # What's the edge called?
    edg.loc <- strsplit(edg.loc, ".", fixed = TRUE)[[1]] # Pull it apart into the node names.

    left.node.idx          <- which(node.names.mat$node.name == edg.loc[1])
    right.node.idx         <- which(node.names.mat$node.name == edg.loc[2])
    edge.mat.loc[i,]       <- c(left.node.idx, right.node.idx)
    edge.mat.names.loc[i,] <- c(edg.loc[1], edg.loc[2])

  }

  # Make adjacency mat from the (node number) edge mat:
  adj.mat.loc <- edges2adj(edge.mat.loc)

  if(standard.potentialsQ == T) {

    # Coefs come out of loglin basically in "Ising" parameterization (cf. Notes). Shift them over to standard parameterization
    crf.obj <- make.empty.field(adj.mat = adj.mat.loc, parameterization.typ="standard", plotQ=T)

    par.count <- 1 # Running index for the parameter (energy) vector
    # Extract and store node parameters and potentials
    for(i in 1:length(node.idxs)) {

      # Remove any decorations:
      a.tau <- as.numeric(model.coefs[node.idxs[i]][[1]]) # an H vector (i.e. node parameters in Ising format, H = (H1, -H1))

      # Put in node parameter (tau) in standard parameterization and then exponentiate to potential
      a.tau                  <- a.tau + a.tau[1] # H is (H1, -H1) so tau_i1 = H1 + H1, tau_i2 = -H1 + H1
      crf.obj$par[par.count] <- a.tau[1]

      a.node.pot           <- exp(a.tau)
      crf.obj$node.pot[i,] <- a.node.pot

      par.count <- par.count + 1
    }

    # Extract and store edge parameters and potentials
    for(i in 1:length(edge.idxs)) {

      # Remove any decorations:
      a.omega <- array(as.numeric(model.coefs[edge.idxs[i]][[1]]), c(2,2)) # a J matrix (i.e. edge parameters in Ising format, J = (J1, -J1), (-J1, J))

      if( round(a.omega[1,2],6) != round(a.omega[2,1],6) ) { # Check: But round first in case they are a little asymmetric
        print(paste("**Edge: ", node.and.edge.names[edge.idxs[i]]))
        print(paste("**Edge# ", i, " potential not symmetric"))
        print(a.omega)
        print(a.omega[1,2])
        print(a.omega[2,1])
        print(a.omega[1,2] - a.omega[2,1])
        stop("** Encountered non-symmetric edge parameter matrix! Cannot write in standard parameterization! CHECK!!")
      }

      # Put omega in standard parameterization and then exponentiate:
      a.omega <- a.omega + a.omega[1,1] # J is ((J1, -J1), (-J1, J1)) so omega_i11 = J1 + J1, omega_i12 = -J1 + J1, omega_i21 = -J1 + J1, omega_i22 = J1 + J1

      crf.obj$par[par.count] <- a.omega[1,1]

      a.edge.pot            <- exp(a.omega)
      crf.obj$edge.pot[[i]] <- a.edge.pot

      par.count <- par.count + 1

    }

  } else {

    # Just put the coefs into the crf object as they came out of loglin() which is "Ising" parameterization (cf. Notes):
    crf.obj <- make.empty.field(adj.mat = adj.mat.loc, parameterization.typ="general", plotQ=plotQ)

    # Node potentials:
    for(i in 1:length(node.idxs)) {

      # Remove any decorations:
      a.tau <- as.numeric(model.coefs[node.idxs[i]][[1]]) # an H vector (i.e. node parameters in Ising format, H = (H1, -H1))

      # Put in node parameter (tau) in as it came out of loglin() parameterization and then exponentiate to potential
      a.node.pot           <- exp(a.tau)
      crf.obj$node.pot[i,] <- a.node.pot
    }

    # Edge potentials:
    for(i in 1:length(edge.idxs)) {

      # Remove any decorations:
      a.omega <- array(as.numeric(model.coefs[edge.idxs[i]][[1]]), c(2,2)) # a J matrix (i.e. edge parameters in Ising format, J = (J1, -J1), (-J1, J))

      # Put in omega as it comes out of loglin() and then exponentiate:
      a.edge.pot            <- exp(a.omega)
      crf.obj$edge.pot[[i]] <- a.edge.pot

    }

    crf.obj$par <- make.par.from.all.potentials(crf.obj)$par

  }
  #dump.crf(crf.obj)

  return(crf.obj)

}


#' XXXX
#'
#' XXXX
#'
#' @details XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
glmnet_logistic_fit2graph_info <- function(a.glmnet.logistic.fit, num.nodes, crf.obj=NULL, lambda="lambda.1se", plotQ=F){

  if(is.null(crf.obj)) { # If no crf.obj container is passed in, use a saturated graph

    if(is.null(num.nodes)) {
      stop("num.nodes in the glmnet fit must be specified if no crf.obj is specified!")
    }

    g.saturated.start           <- erdos.renyi.game(num.nodes, 1, typ="gnp")
    adj.saturated.loc           <- as.matrix(as_adj(g.saturated.start))
    colnames(adj.saturated.loc) <- 1:num.nodes
    rownames(adj.saturated.loc) <- 1:num.nodes
    a.crf.loc                   <- make.empty.field(adj.mat = adj.saturated.loc, parameterization.typ = "standard", plotQ = plotQ)

  } else {
    # ****NOTE: This assumes the glmnet fit was done for a model with a graph implied by the crf.obj graph!
    a.crf.loc <- crf.obj
  }

  num.nodes.loc <- a.crf.loc$n.nodes

  loc.fit.coefs <- as.numeric(coef(a.glmnet.logistic.fit, s = lambda))

  # Indices of parameters with non-zero values. We use these below to get the non-zero edge parameters
  # **NOTE: Assumes no intercept fit so subtract 1 from non-zero indices
  nonz.idxs <- which(loc.fit.coefs != 0 ) - 1

  # Find which nodes/edges are associated with which (non-zero) params. **NOTE, assumes parameterization is standard pairwise
  prms2nds.lst <- params2nodes.list(a.crf.loc)[nonz.idxs]
  print(prms2nds.lst)

  prm2edge.info.mat <- NULL
  for(i in 1:length(prms2nds.lst)) {

    if(length(prms2nds.lst[[i]]) == 2) { # Grab parameters for the edges. These are the important ones for the model structure (i.e. the glmnet selected graph).
      prm2edge.vec <- c(nonz.idxs[i], loc.fit.coefs[ nonz.idxs[i]+1 ], prms2nds.lst[[i]])
      #print(prm2edge.vec)

      prm2edge.info.mat <- rbind(prm2edge.info.mat, prm2edge.vec)
    }

  }
  colnames(prm2edge.info.mat) <- c("param", "edge.coef", "left.edge.node", "right.edge.node")
  rownames(prm2edge.info.mat) <- NULL
  #print(prm2edge.info.mat)

  # Note: glmnet may have forced some of the node parameters to be 0, but we always want to keep all node parameters.
  # In standard parameterization, these are the parameters associated with the nodes:
  node.param.idxs <- 1:num.nodes.loc
  prm2node.info.mat <- cbind(node.param.idxs, loc.fit.coefs[ node.param.idxs+1 ], node.param.idxs)
  colnames(prm2node.info.mat) <- c("param", "node.coef", "node")
  #print(prm2node.info.mat)

  if(plotQ==T) {
    # Recover the graph:
    adj.recovered.loc <- edges2adj(prm2edge.info.mat[,c(3,4)], n.nodes = num.nodes.loc)
    print(adj.recovered.loc)

    #grph.recovered.loc <- adj2formula(adj.recovered.loc, edge.only.formulaQ = T, Xoption = F)
    #gp.recovered.loc   <- ug(grph.recovered.loc, result = "igraph")
    gp.recovered.loc <- graph_from_adjacency_matrix(adj.recovered.loc, mode = "undirected")

    # if(length(dev.list()) !=0) {
    #   dev.off()
    # }
    plot(gp.recovered.loc)
  }

  par.vec <- c(prm2node.info.mat[,2], prm2edge.info.mat[,2])

  param.info.list <- list(prm2node.info.mat, prm2edge.info.mat, par.vec)
  names(param.info.list) <- c("node.parameters.mat", "edge.parameters.mat", "parameter.vector")

  return(param.info.list)

}
