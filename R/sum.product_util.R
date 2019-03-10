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
get.root.paths <- function(graph.eq, root.node) {

  und.gph    <- ug(graph.eq)
  leaf.nodes <- leaves(und.gph)

  all.root.pths <- lapply(1:length(leaf.nodes), function(xx){get.path(und.gph, root.node, leaf.nodes[xx])})

  # Check and see if any of the paths are NULL. Probably chose a leaf node as root in that case:
  null.pthsQ    <- sapply(1:length(all.root.pths), function(xx){is.null(all.root.pths[[xx]])})
  drp.pths <- which(null.pthsQ == TRUE)
  if(length(drp.pths) !=0) {
    all.root.pths <- all.root.pths[-drp.pths]
  }

  return(all.root.pths)

}


#' Process a list of paths into serial schedules (forward and backward)
#'
#' Process a list of paths into serial schedules (forward and backward)
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
paths.to.serial.scheds <- function(paths.list) {

  num.paths  <- length(paths.list)
  path.lengs <- sapply(1:num.paths, function(xx){length(paths.list[[xx]])})

}
