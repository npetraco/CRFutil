#' Plot a a sample of states node-wise
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
