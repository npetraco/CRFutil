#' marginal edge empirical probs
#'
#' Empirical Pr(X1), Pr(X2), Pr(X1,X2), Pr(X1|X2), Pr(X2|X1) estimates from edge data
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bayes.bels.plot <- function(posterior.edge.belief.info, type="X1|X2", ymax=2500, edge.empirical.prob.info=NULL){


  distc.nmes <- names(posterior.edge.belief.info)[c(4,5)] # Edge conditional distributions
  distm.nmes <- names(posterior.edge.belief.info)[c(2,3)] # Edge marginal node distributions

  #print(distc.nmes)
  #print(distm.nmes)

  # Choose conditional direction:
  if(type=="X1|X2") {
    distc.nme <- distc.nmes[1]
    distm.nme <- distm.nmes[1]
  } else if(type=="X2|X1") {
    distc.nme <- distc.nmes[2]
    distm.nme <- distm.nmes[2]
  } else {
    stop("type must = X1|X2 or X2|X1")
  }

  cond.sect.nmes <- colnames(posterior.edge.belief.info[[distc.nme]])
  marg.sect.nmes <- colnames(posterior.edge.belief.info[[distm.nme]])
  print(distc.nme)
  print(cond.sect.nmes)
  print(distm.nme)
  print(marg.sect.nmes)

  # Shut off plot window if open
  if(!is.null(dev.list())){
    dev.off()
  }
  par(mfrow=c(2,2)) # split into 4

  # sector 1,1
  cond.sect.nme <- cond.sect.nmes[1]
  marg.sect.nme <- marg.sect.nmes[1]
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax))
  abline(v=med.lne, lwd=4)
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 1,2
  cond.sect.nme <- cond.sect.nmes[2]
  marg.sect.nme <- marg.sect.nmes[1]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax))
  abline(v=med.lne, lwd=4)
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 2,1
  cond.sect.nme <- cond.sect.nmes[3]
  marg.sect.nme <- marg.sect.nmes[2]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax))
  abline(v=med.lne, lwd=4)
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 2,2
  cond.sect.nme <- cond.sect.nmes[4]
  marg.sect.nme <- marg.sect.nmes[2]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax))
  abline(v=med.lne, lwd=4)
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  print(paste0("If X1_||_X2 then: ",distc.nmes[1], " = ", distm.nmes[1]))
  print(paste0("                  ",distc.nmes[2], " = ", distm.nmes[2]))

}
