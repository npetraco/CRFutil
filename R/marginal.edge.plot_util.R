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
marginal.edge.bayes.bels.plot <- function(posterior.edge.belief.info, type="X1|X2", ymax=2500, prob.level=0.95){


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
  #print(distc.nme)
  #print(cond.sect.nmes)
  #print(distm.nme)
  #print(marg.sect.nmes)

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
  cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
  marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)

  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax), col="red")
  abline(v=med.lne, lwd=4)
  points(cond.pi, c(0,0), pch=8, col="firebrick4")
  points(marg.pi, c(0,0), pch=16, col="black")
  med.in.intQ <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
  if(med.in.intQ == TRUE) {
    t.nme <- paste0(t.nme,"**")
  }
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 1,2
  cond.sect.nme <- cond.sect.nmes[2]
  marg.sect.nme <- marg.sect.nmes[1]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
  marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)

  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax), col="red")
  abline(v=med.lne, lwd=4)
  points(cond.pi, c(0,0), pch=8, col="firebrick4")
  points(marg.pi, c(0,0), pch=16, col="black")
  med.in.intQ <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
  if(med.in.intQ == TRUE) {
    t.nme <- paste0(t.nme,"**")
  }
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 2,1
  cond.sect.nme <- cond.sect.nmes[3]
  marg.sect.nme <- marg.sect.nmes[2]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
  marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)

  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax), col="red")
  abline(v=med.lne, lwd=4)
  points(cond.pi, c(0,0), pch=8, col="firebrick4")
  points(marg.pi, c(0,0), pch=16, col="black")
  med.in.intQ <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
  if(med.in.intQ == TRUE) {
    t.nme <- paste0(t.nme,"**")
  }
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  # sector 2,2
  cond.sect.nme <- cond.sect.nmes[4]
  marg.sect.nme <- marg.sect.nmes[2]
  plt.main <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  plt.main      <- paste0("Bel(",cond.sect.nme,") vs. Med[Bel(",marg.sect.nme,")]")
  x.nme         <- paste0("Bel(",cond.sect.nme,")")
  t.nme         <- paste0("Med[Bel(",marg.sect.nme,")]")
  med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
  cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
  marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)

  hist(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], xlab=x.nme, main=plt.main, xlim=c(0,1), ylim=c(0,ymax), col="red")
  abline(v=med.lne, lwd=4)
  points(cond.pi, c(0,0), pch=8, col="firebrick4")
  points(marg.pi, c(0,0), pch=16, col="black")
  med.in.intQ <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
  if(med.in.intQ == TRUE) {
    t.nme <- paste0(t.nme,"**")
  }
  text(label=t.nme, x=med.lne, y=ymax, adj=1.2)

  print(paste0("If X1_||_X2 then: ",distc.nmes[1], " = ", distm.nmes[1]))
  print(paste0("                  ",distc.nmes[2], " = ", distm.nmes[2]))

}

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
marginal.edge.bayes.plot <- function(posterior.edge.pot.mat, type="pot", prob.level=0.95){

  if(type == "pot") { # for poentials
    tau1.loc    <- posterior.edge.pot.mat[,1]
    tau2.loc    <- posterior.edge.pot.mat[,2]
    omeg.loc    <- posterior.edge.pot.mat[,3]
    xlabs       <- c(expression(e^{tau[1]}), expression(e^{tau[2]}), expression(e^{omega}))
    mains       <- c("exp(tau1)", "exp(tau2)", "exp(omega)")
    indep.indic <- 1
  } else if(type == "logpot") { # for energies
    tau1.loc    <- log(posterior.edge.pot.mat[,1])
    tau2.loc    <- log(posterior.edge.pot.mat[,2])
    omeg.loc    <- log(posterior.edge.pot.mat[,3])
    xlabs       <- c(expression(tau[1]), expression(tau[2]), expression(omega))
    mains       <- c("tau1", "tau2", "omega")
    indep.indic <- 0
  } else {
    stop("Choose: pot or logpot for type!")
  }

  tau1.int   <- HPDI(samples = tau1.loc, prob = prob.level)
  tau2.int   <- HPDI(samples = tau2.loc, prob = prob.level)
  omeg.int   <- HPDI(samples = omeg.loc, prob = prob.level)


  # Shut off plot window if open
  if(!is.null(dev.list())){
    dev.off()
  }
  par(mfrow=c(2,2)) # split into 4

  hist(tau1.loc, xlab=xlabs[1], main=mains[1])
  points(tau1.int, c(0,0), pch=16, col="blue")
  abline(v=indep.indic, lwd=4, col="green")

  hist(tau2.loc, xlab=xlabs[2], main=mains[2])
  abline(v=indep.indic, lwd=4, col="green")
  points(tau2.int, c(0,0), pch=16, col="blue")

  hist(omeg.loc, xlab=xlabs[3], main=mains[3])
  abline(v=indep.indic, lwd=4, col="green")
  points(omeg.int, c(0,0), pch=16, col="blue")

}
