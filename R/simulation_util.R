#' Simulate an MRF with random unary (node) and pair (edge) potentials:
#'
#' Simulate an MRF with random unary (node) and pair (edge) potentials:
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
sim.field.random <- function(adjacentcy.matrix, num.states, num.sims, seed=NULL) {

  mrf.sim.model <- make.crf(adjacentcy.matrix, num.states)

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

  mrf.model.samples <- sample.junction(mrf.sim.model, num.sims)
  colnames(mrf.model.samples) <- as.character(1:num.nodes)
  mrf.sim.info <- list(
    mrf.sim.model,
    mrf.model.samples
  )
  names(mrf.sim.info) <- c("model", "samples")

  return(mrf.sim.info)

}
