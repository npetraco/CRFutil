#' Simulate an MRF with random unary (node) and pair (edge) potentials using a specified adjacency matrix.
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials using a specified adjacency matrix.
#'
#' Uses standard parameterization w11=w22, w12=w21
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
sim.field.random <- function(adjacentcy.matrix, num.states, num.sims, seed=NULL) {

  mrf.sim.model <- make.crf(adjacentcy.matrix, num.states)

  # May be needed for other functionality, so include explicit (standard) parameterization:
  mrf.sim.model <- make.features(mrf.sim.model)
  mrf.sim.model <- make.par(mrf.sim.model, nrow(mrf.sim.model$node.pot) + length(mrf.sim.model$edge.pot))
  for(i in 1:nrow(mrf.sim.model$node.par)){
    mrf.sim.model$node.par[i,1,1] <- i
  }
  for(i in 1:length(mrf.sim.model$edge.par)){
    mrf.sim.model$edge.par[[i]][1,1,1] <- nrow(mrf.sim.model$node.pot) + i
    mrf.sim.model$edge.par[[i]][2,2,1] <- nrow(mrf.sim.model$node.pot) + i
  }

  num.nodes <- nrow(adjacentcy.matrix)

  # Make up random node and edge potentials:
  if(!is.null(seed)){
    set.seed(seed)
  }
  rand.pots <- runif(num.nodes + mrf.sim.model$n.edges)

  # Node potentials
  pos <- rand.pots[1:num.nodes]
  neg <- 1-pos

  mrf.sim.model$node.pot <- cbind(pos,neg)

  # Edge potentials
  edge.pots.vec <- rand.pots[(num.nodes+1):length(rand.pots)]
  for (i in 1:mrf.sim.model$n.edges) {
    # Make symmetric edge potentials
    trans.1122 <- edge.pots.vec[i]
    trans.prob <- rbind(
      c(trans.1122, 1-trans.1122),
      c(1-trans.1122, trans.1122)
    )

    mrf.sim.model$edge.pot[[i]] <- trans.prob
  }

  # Extract parameter vector:
  mrf.sim.model$par <- make.par.from.potentials(mrf.sim.model)

  # Scale the potentials to conform with the parameter vector:
  rescaled.pots     <- make.pots(mrf.sim.model$par, mrf.sim.model, rescaleQ=FALSE, replaceQ=TRUE, printQ=FALSE)

  #What nodes are associated with what parameter?
  nodes2parameters  <- nodes2params.list(mrf.sim.model, storeQ = TRUE)


  mrf.model.samples <- sample.junction(mrf.sim.model, num.sims)
  #colnames(mrf.model.samples) <- as.character(1:num.nodes)
  colnames(mrf.model.samples) <- as.character(colnames(adjacentcy.matrix)) # How do we know for sure these are the column names??
  mrf.sim.info <- list(
    mrf.sim.model,
    mrf.model.samples
  )
  names(mrf.sim.info) <- c("model", "samples")

  return(mrf.sim.info)

}

#' Simulate an MRF with random unary (node) and pair (edge) potentials using a specified adjacency matrix.
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials using a specified adjacency matrix.
#'
#' Uses simplest standard parameterization I can think of: tau1 = tau2 ..., w11=w22, w12=w21, ...
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
sim.field.random.simplest <- function(adjacentcy.matrix, num.states, num.sims, seed=NULL) {

  mrf.sim.model <- make.crf(adjacentcy.matrix, num.states)

  # May be needed for other functionality, so include explicit (standard) parameterization:
  mrf.sim.model <- make.features(mrf.sim.model)
  mrf.sim.model <- make.par(mrf.sim.model, 2)
  for(i in 1:nrow(mrf.sim.model$node.par)){
    mrf.sim.model$node.par[i,1,1] <- 1
  }
  for(i in 1:length(mrf.sim.model$edge.par)){
    #mrf.sim.model$edge.par[[i]][1,1,1] <- nrow(mrf.sim.model$node.pot) + i
    #mrf.sim.model$edge.par[[i]][2,2,1] <- nrow(mrf.sim.model$node.pot) + i
    mrf.sim.model$edge.par[[i]][1,1,1] <- 2
    mrf.sim.model$edge.par[[i]][2,2,1] <- 2
  }


  # Make up random node weights:
  num.nodes <- nrow(adjacentcy.matrix)
  if(!is.null(seed)){
    set.seed(seed)
  }
  pos <- runif(num.nodes)
  neg <- 1-pos

  # ADD facility for non random node weights

  mrf.sim.model$node.pot <- cbind(pos,neg)

  for (i in 1:mrf.sim.model$n.edges) {
    # Make up random symmetric edge weights
    if(!is.null(seed)){
      set.seed(seed)
    }
    trans.1122 <- runif(1)
    trans.prob <- rbind(
      c(trans.1122, 1-trans.1122),
      c(1-trans.1122, trans.1122)
    )

    mrf.sim.model$edge.pot[[i]] <- trans.prob
  }

  # Extract parameter vector:
  mrf.sim.model$par <- make.par.from.potentials(mrf.sim.model)

  # Scale the potentials to conform with the parameter vector:
  rescaled.pots     <- make.pots(mrf.sim.model$par, mrf.sim.model, rescaleQ=FALSE, replaceQ=TRUE, printQ=FALSE)

  #What nodes are associated with what parameter?
  nodes2parameters  <- nodes2params.list(mrf.sim.model, storeQ = TRUE)


  mrf.model.samples <- sample.junction(mrf.sim.model, num.sims)
  colnames(mrf.model.samples) <- as.character(1:num.nodes)
  mrf.sim.info <- list(
    mrf.sim.model,
    mrf.model.samples
  )
  names(mrf.sim.info) <- c("model", "samples")

  return(mrf.sim.info)

}


#' Simulate an random MRF with random unary (node) and pair (edge) potentials.
#'
#' Simulate an random MRF with random unary (node) and pair (edge) potentials.
#'
#' Uses standard parameterization w11=w22, w12=w21. Graph is generated with erdos.renyi.game
#' in igraph. If an edge appears in the graph, a relatively high to high potential (i.e. non-zero)
#' should be generated to ensure the edge is potentially strong enough to be found from a random
#' sample from the graph. However, the user can choose the the general strength and variability of
#' the edge and node potentials though the edge.pot.mean/sd and node.pot.mean/sd arguments. They
#' are generated with a t-distribution. The default nu is 30 making the random selections approximately
#' normal. The number of edges in a randomly generated graph can be chosen by specification or by
#' probability.
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
sim.random.field <- function(adj.mat=NULL, num.nodes, num.edges=NULL, prob.edge=1, node.pot.nu=30, node.pot.mean, node.pot.sd, edge.pot.nu=30, edge.pot.mean, edge.pot.sd, num.sims, seed.node.pots=NULL, seed.edge.pots=NULL, seed.sample=NULL, plotQ=F) {

  if(is.null(adj.mat)) { # If no adjacency matrix is specified, generate one

    num.nodes.loc <- num.nodes

    if(!is.null(num.edges)) {
      graph.loc <- erdos.renyi.game(n = num.nodes.loc, p.or.m = num.edges, type = "gnm", directed = F, loops = F)
    } else {
      graph.loc <- erdos.renyi.game(n = num.nodes.loc, p.or.m = prob.edge, type = "gnp", directed = F, loops = F)
    }

    adj.mat.loc           <- as.matrix(as_adj(graph.loc))
    colnames(adj.mat.loc) <- 1:num.nodes.loc
    rownames(adj.mat.loc) <- 1:num.nodes.loc

  } else {                  # If an adjacency matric is specified, use that.
    num.nodes.loc <- num.nodes
    adj.mat.loc   <- adj.mat
    colnames(adj.mat.loc) <- 1:num.nodes.loc # Id adjacency matrix has row/col, names, just number them.
    rownames(adj.mat.loc) <- 1:num.nodes.loc

    graph.loc <- graph_from_adjacency_matrix(adjmatrix = adj.mat.loc, mode = "undirected")
  }

  # if(plotQ==T) {
  #   plot(graph.loc)
  # }

  model.loc <- make.empty.field(adj.mat = adj.mat.loc, parameterization.typ = "standard", plotQ = plotQ)
  #print(model.loc$n.nodes)
  #print(model.loc$n.edges)

  if(!is.null(seed.node.pots)){
    set.seed(seed.node.pots)
  }
  # Node pots/pars. Abs the pots to make sure they're positive
  node.pot.vec.loc <- abs(ggdist::rstudent_t(n = model.loc$n.nodes, df = node.pot.nu, mu = node.pot.mean, sigma = node.pot.sd))
  node.par.loc     <- log(node.pot.vec.loc)

  if(!is.null(seed.edge.pots)){
    set.seed(seed.edge.pots)
  }
  # Edge pots/pars. Abs the pots to make sure they're positive
  edge.pot.vec.loc <- abs(ggdist::rstudent_t(n = model.loc$n.edges, df = edge.pot.nu, mu = edge.pot.mean, sigma = edge.pot.sd))
  edge.par.loc     <- log(edge.pot.vec.loc)

  model.loc$par <- c(node.par.loc, edge.par.loc)
  out.pots.loc  <- make.pots(parms = model.loc$par, crf = model.loc, replaceQ = T)
  #dump.crf(model.loc)

  if(!is.null(seed.sample)){
    set.seed(seed.sample)
  }
  model.loc.samples <- sample.junction(model.loc, num.sims)
  colnames(model.loc.samples) <- as.character(1:num.nodes.loc)
  sim.info <- list(
    model.loc,
    model.loc.samples
  )
  names(sim.info) <- c("model", "samples")

  return(sim.info)

}


#' Simulate an MRF with random unary (node) and pair (edge) potentials.
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials.
#'
#' Uses xx standard parameterization
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
sim.ising.field.random <- function(adjacentcy.matrix, parameterization.type="a", num.sims, seed=NULL) {

  mrf.sim.model <- make.crf(adjacentcy.matrix, num.states)

  # Make up random node weights:
  num.nodes <- nrow(adjacentcy.matrix)
  if(!is.null(seed)){
    set.seed(seed)
  }
  pos <- runif(num.nodes)
  neg <- 1-pos
  mrf.sim.model$node.pot <- cbind(pos,neg)


  for (i in 1:mrf.sim.model$n.edges) {
    # Make up random symmetric edge weights
    if(!is.null(seed)){
      set.seed(seed)
    }

    # Standard parameterization a. w11=w22, w12=w21=1-w11
    if(parameterization.type=="a") {
      trans.1122 <- runif(1)
      trans.prob <- rbind(
        c(trans.1122, 1-trans.1122),
        c(1-trans.1122, trans.1122)
      )
    }

    # Standard parameterization b. w11 != w22, w12=w21
    if(parameterization.type=="b") {
      pots <- runif(3)
      trans.prob <- rbind(
        c(pots[1], pots[3]),
        c(pots[3], pots[2])
      )
    }

    # Standard parameterization c. w11 = 1-w12, w22=1-w21
    if(parameterization.type=="c") {
      pots <- runif(2)
      trans.prob <- rbind(
        c(  pots[1], 1-pots[2]),
        c(1-pots[1], pots[2]  )
      )
    }

    # Standard parameterization d. w11 != w12 != w22 != w21
    if(parameterization.type=="d") {
      trans.prob <- array(runif(4), c(2,2))
    }

    mrf.sim.model$edge.pot[[i]] <- trans.prob
  }

  mrf.model.samples <- sample.junction(mrf.sim.model, num.sims)
  colnames(mrf.model.samples) <- as.character(1:num.nodes)
  mrf.sim.info <- list(
    mrf.sim.model,
    mrf.model.samples
  )
  names(mrf.sim.info) <- c("model", "samples")

  return(mrf.sim.info)

}
