#' Empirical joint distribution
#'
#' Empirical joint distribution
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
fit_empirical <- function(samples) {

  # Take a look at the empirical frequency-counts of the states:
  configs.and.counts <- as.data.frame(ftable(data.frame(samples)))
  freq.idx           <- ncol(configs.and.counts)
  config.counts      <- configs.and.counts[,freq.idx]

  # Estimate configuration state probabilities is by the EMPIRICAL relative frquencies:
  config.freqs <- config.counts/sum(config.counts)

  # Order nodes in increasing order (in case they got out of wack):
  # NOTE: these should only be numbers for testing purposes!
  # BUT data.frame puts an X in front. Remove it to order
  node.nums             <- unlist(strsplit(colnames(configs.and.counts)[-freq.idx], split = "X"))
  node.nums             <- node.nums[-which(node.nums == "")]
  node.col.reorder.idxs <- order(node.nums)

  #configs.and.counts           <- configs.and.counts[, node.col.reorder.idxs]
  configs.and.counts                     <- cbind(configs.and.counts[-freq.idx],config.freqs)
  colnames(configs.and.counts)[freq.idx] <- "Emp.Freqs"

  return(configs.and.counts)

}
