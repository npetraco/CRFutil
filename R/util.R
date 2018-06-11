#' Log sum exp trick. From Brendon Brewer's DNest code:
#'
#' Handy for calculating Z from a vector of log potentials.
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
logsumexp <- function(logv) {
  n <- length(logv)
  max.logv <- max(logv)

  answer <- 0

  for(i in 1:n){
    answer <- answer + exp(logv[i] - max.logv)
  }
  answer <- max.logv + log(answer);

  return(answer)

}


#' Log sum exp trick. From Brendon Brewer's DNest code:
#'
#' Handy for calculating Z from a vector of log potentials.
#' A less readable but shorter Log sum exp:
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
logsumexp2 <- function(logv)
{
  n <- length(logv)
  max.logv <- max(logv)

  answer <-  max.logv + log(cumsum(c(0,exp(logv - max.logv)))[n+1])

  return(answer)

}

#' Code from prodlim library to match a row in a matrix
#'
#' XXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
row.match <- function (x, table, nomatch = NA)
{
  if (class(table) == "matrix")
    table <- as.data.frame(table)
  if (is.null(dim(x)))
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}

