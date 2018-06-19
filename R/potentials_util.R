#' Shift node potentials back to the originally fit parameter vector
#'
#' For testing purposes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
shift.pots <- function(crf) {

  parms <- crf$par

  # Loop over elements of parameter index matrix (also called nodeMap and edgeMap in UGM)
  # and put elements of exp(parms) where they belong

  node.par         <- crf$node.par[,,1]
  node.pot.shifted <- array(1, c(dim(node.par),1) )
  for(i in 1:nrow(node.par)) {
    for(j in 1:ncol(node.par)) {
      if(node.par[i,j] != 0) {
        par.idx <- node.par[i,j]
        node.pot.shifted[i,j,1] <- exp(parms[par.idx])
      }
    }
  }


  edge.par <- crf$edge.par
  edge.pot.shifted <- rep(list(array(1,c(2,2))), length(edge.par))
  for(k in 1:length(edge.par)) {
    for(i in 1:nrow(edge.par[[k]])) {
      for(j in 1:ncol(edge.par[[k]])) {
        if(edge.par[[k]][i,j,1] != 0) {
          par.idx <- edge.par[[k]][i,j,1]
          edge.pot.shifted[[k]][i,j] <- exp(parms[par.idx])
        }
      }
    }
  }

  crf$node.pot <- node.pot.shifted
  crf$edge.pot <- edge.pot.shifted
  print("Potentials shifted to parameter vector via node.par and edge.par")

}


#' Make stand alone potentials from input parameter vector.
#'
#' Does same thing as shift.pots but returns potentials instead of putting them into the CRF object.
#' Eventually this could be combined with shift.pots
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.pots <- function(parms, crf) {

  # Loop over elements of parameter index matrix (also called nodeMap and edgeMap in UGM)
  # and put elements of exp(parms) where they belong
  node.par         <- crf$node.par[,,1]
  node.pot.shifted <- array(1, c(dim(node.par),1) )
  for(i in 1:nrow(node.par)) {
    for(j in 1:ncol(node.par)) {
      if(node.par[i,j] != 0) {
        par.idx <- node.par[i,j]
        node.pot.shifted[i,j,1] <- exp(parms[par.idx])
      }
    }
  }

  edge.par <- crf$edge.par
  edge.pot.shifted <- rep(list(array(1,c(2,2))), length(edge.par))
  for(k in 1:length(edge.par)) {
    for(i in 1:nrow(edge.par[[k]])) {
      for(j in 1:ncol(edge.par[[k]])) {
        if(edge.par[[k]][i,j,1] != 0) {
          par.idx <- edge.par[[k]][i,j,1]
          edge.pot.shifted[[k]][i,j] <- exp(parms[par.idx])
        }
      }
    }
  }

  # Don't worry about the .shifted extension. It's just so I could copy code from shift.pots.
  return(list(node.pot.shifted,edge.pot.shifted))


}
