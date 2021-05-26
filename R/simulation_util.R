#' Simulate an MRF with random unary (node) and pair (edge) potentials.
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials.
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

#' Simulate an MRF with random unary (node) and pair (edge) potentials.
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials.
#'
#' Uses simplest standard parameterization I can tink of: tau1 = tau2 ..., w11=w22, w12=w21, ...
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
