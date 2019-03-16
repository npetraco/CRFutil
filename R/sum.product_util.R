#' Get a (shortest) path between two nodes
#'
#' Get a (shortest) path between two nodes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.path <- function(a.graphnel.obj, start.node, end.node){

  path.info <- sp.between(a.graphnel.obj, as.character(start.node), as.character(end.node))[[1]]
  if(is.na(path.info[1])) {
    #print("No path!")
    sepath <- NULL
  } else {
    sepath <- path.info$path_detail
  }

  return(sepath)
}


#' Get all paths from chosen root node to all leaf nodes
#'
#' Get all paths from chosen root node to all leaf nodes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.root.paths <- function(graph.obj, root.node, serial.schedsQ=FALSE) {

  if(class(graph.obj) == "formula") {
    und.gph    <- ug(graph.obj)
  } else if(class(graph.obj) == "graphNEL") {
    und.gph    <- graph.obj
  } else {
    stop("Must enter a formula or graphNEL object for the graph.obj arguement.")
  }
  leaf.nodes <- leaves(und.gph)

  all.root.pths <- lapply(1:length(leaf.nodes), function(xx){get.path(und.gph, root.node, leaf.nodes[xx])})

  # Check and see if any of the paths are NULL. Probably chose a leaf node as root in that case:
  null.pthsQ <- sapply(1:length(all.root.pths), function(xx){is.null(all.root.pths[[xx]])})
  drp.pths   <- which(null.pthsQ == TRUE)
  if(length(drp.pths) !=0) {
    all.root.pths <- all.root.pths[-drp.pths]
  }

  if(serial.schedsQ == TRUE) {
    sched.mat.forward  <- array(NA,c(length(all.root.pths), max(sapply(1:length(all.root.pths), function(xx){length(all.root.pths[[xx]])}))))
    sched.mat.backward <- array(NA,c(length(all.root.pths), max(sapply(1:length(all.root.pths), function(xx){length(all.root.pths[[xx]])}))))
    for(i in 1:length(all.root.pths)) {
      rev.pth <- rev(all.root.pths[[i]])
      pth     <- all.root.pths[[i]]
      for(j in 1:length(all.root.pths[[i]])) {
        sched.mat.forward[i,j]  <- rev.pth[j]
        sched.mat.backward[i,j] <- pth[j]
      }
    }
    all.root.pths <- list(sched.mat.forward, sched.mat.backward)
    names(all.root.pths) <- c("forwrd", "backward")
  }

  return(all.root.pths)

}

#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf2pwfg <- function(graph.obj, plotQ=FALSE) {

  if(class(graph.obj) == "formula") {
    und.gph    <- ug(graph.obj)
  } else if(class(graph.obj) == "graphNEL") {
    und.gph    <- graph.obj
  } else {
    stop("Must enter a formula or graphNEL object for the graph.obj arguement.")
  }

  # Use igraph functionality to build pair-wise factor graph from UG model
  adjmat <- as(und.gph, "matrix")
  adjmat.uptri <- upper.tri(adjmat) * adjmat

  # edges:
  edgmat <- which(adjmat.uptri == 1, arr.ind = T)
  edgmat <- cbind(
    as.numeric(rownames(adjmat)[edgmat[,1]]),
    as.numeric(colnames(adjmat)[edgmat[,2]])
  )
  edgmat <- t(sapply(1:nrow(edgmat), function(xx){sort(edgmat[xx,])}))

  edgmat2 <- NULL
  # add egde factors between connected nodes:
  for(i in 1:nrow(edgmat)) {
    factor.nme <- paste0("f",edgmat[i,1],edgmat[i,2])
    edgmat2 <- rbind(
      edgmat2,
      c(edgmat[i,1], factor.nme),
      c(factor.nme, edgmat[i,2])
    )
  }

  # tack on node factors:
  edgmat2 <- rbind(
    cbind(
      paste0("f",sort(und.gph@nodes)),
      sort(und.gph@nodes)
    ),
    edgmat2
  )

  fg <- graph_from_data_frame(data.frame(edgmat2), directed = FALSE)
  nde.nms <- V(fg)$name
  V(fg)$type <- sapply(1:length(nde.nms), function(xx){length(strsplit(nde.nms[xx], split = "f")[[1]])})

  if(plotQ==TRUE){
    # Plot factor graph:
    cols <- c("steelblue", "red")
    shps <- c("circle", "square")
    plot(fg,
         vertex.color = cols[as.numeric(V(fg)$type)],
         vertex.shape = shps[as.numeric(V(fg)$type)]
    )
  }

  # Convert back to graphNEL format:
  fg <- igraph.to.graphNEL(fg)
  return(fg)
}


#' Make factor to variable (node) message
#'
#' Make factor to variable (node) message
#' \mu_{f\rightarrow X} = \sum_{\backslash X} \Big( f(X) \prod_{Y \in \text{ne}(f)\backslash X} \mu_{Y\rightarrow f} \Big)
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.f2v.msg <- function(in.v.msgs.list = NULL, f.msg, out.v.nme){

  if(is.null(in.v.msgs.list)) {
    # v.nme check for initalization??
    msg.out <- f.msg # For initialization case
  } else {
    msg.prod <- f.msg %a*% ar_prod_list(in.v.msgs.list)
    msg.out <- ar_marg(msg.prod, out.v.nme)
  }

  return(msg.out)

}


#' Make factor to variable (node) message
#'
#' Make variable (node) to factormessage
#' XXXXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.v2f.msg <- function(in.f.msgs.list = NULL){

  if(is.null(in.f.msgs.list)) {
    msg.out <- "id" # For initialization case
  } else {
    msg.out <- ar_prod_list(in.f.msgs.list)
  }

  return(msg.out)


}
