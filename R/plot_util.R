#' Plot a sample of states node-wise
#'
#' XX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf.sample.plot <- function(samples) {

  num.nodes <- ncol(samples)
  num.samps <- nrow(samples)

tabl <- table(as.numeric(samples), as.numeric(gl(num.nodes, num.samps)), dnn=c("states", "nodes"))
print(addmargins(tabl))

barplot(tabl, ylab="Frequency", xlab="nodes", main="State Distributions", col=c("green", "red" ), beside=TRUE, width=.3)
#legend("right", title="States", legend= c(1,2), fill =c("turquoise4", "turquoise2"), box.lty=0)

}


#' Plot a crf object
#'
#' XX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
plot_crf <- function(crf.obj, type="g", samples=NULL) {

  adj.loc <- edges2adj(edge.mat = crf.obj$edges, n.nodes = crf.obj$n.nodes, plotQ=FALSE)
  #print(adj.loc)

  grph.loc <- graph_from_adjacency_matrix(adj.loc, mode = "undirected")

  if(type=="g"){

    plot(grph.loc)

  } else if(type=="ge") {

    pars.loc      <- crf.obj$par
    node.pars.loc <- pars.loc[1:crf.obj$n.nodes]
    edge.pars.loc <- pars.loc[(crf.obj$n.nodes + 1):length(pars.loc)]

    par(mar = c(2, 2, 2, 2))
    par(mfrow=c(2,2))
    plot(node.pars.loc, typ="h", main="node pars", ylab="", xlab="", ylim=c(min(node.pars.loc), max(node.pars.loc)) )
    plot(edge.pars.loc, typ="h", main="edge pars", ylab="", xlab="", ylim=c(min(edge.pars.loc), max(edge.pars.loc)) )

    plot_crf(crf.obj)

    if(!is.null(samples)){
      mrf.sample.plot(samples)
    }

    edge.mat.loc           <- cbind(1:crf.obj$n.edges, crf.obj$edges)
    colnames(edge.mat.loc) <- c("edge.num", "left.node", "right.node")
    print(edge.mat.loc)

    # Reset screen parameters
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    par(mfrow=c(1,1))

  } else if(type=="gp") {

    pots.loc      <- exp(crf.obj$par)
    node.pots.loc <- pots.loc[1:crf.obj$n.nodes]
    edge.pots.loc <- pots.loc[(crf.obj$n.nodes + 1):length(pots.loc)]

    par(mar = c(2, 2, 2, 2))
    par(mfrow=c(2,2))
    plot(node.pots.loc, typ="h", main="node pots", ylab="", xlab="", ylim=c(0, max(node.pots.loc)) )
    plot(edge.pots.loc, typ="h", main="edge pots", ylab="", xlab="", ylim=c(0, max(edge.pots.loc)) )

    plot_crf(crf.obj)

    if(!is.null(samples)){
      mrf.sample.plot(samples)
    }

    edge.mat.loc           <- cbind(1:crf.obj$n.edges, crf.obj$edges)
    colnames(edge.mat.loc) <- c("edge.num", "left.node", "right.node")
    print(edge.mat.loc)

    # Reset screen parameters
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    par(mfrow=c(1,1))

  } else {
    stop("type should be g (graph only), ge (graph with parameters) or gp (graph with potentials)")
  }

}
