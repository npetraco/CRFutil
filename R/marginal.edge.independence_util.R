#' marginal edge bayes independence indicators
#'
#' Outputs statistics on posteriors for tau1, tau2, omega
#' Handy for checking for independencies between X1 and X2
#'
#' @param rescaled.posterior.potentials Rescaled posterior potentials from marginal.edge.bayes.loglin.
#'
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bayes.potentials.indep.indic <- function(rescaled.posterior.potentials, node.names=NULL, state.names=NULL, prob.level=0.95, printQ=FALSE){

  if(is.null(node.names)) {
    loc.node.names <- c("X1", "X2")
  } else {
    loc.node.names <- node.names
  }

  if(is.null(state.names)) {
    loc.state.names <- c("1","2")
  } else {
    loc.state.names <- state.names
  }

  # Posterior potentials
  pot.tau1.loc <- rescaled.posterior.potentials[,1]
  pot.tau2.loc <- rescaled.posterior.potentials[,2]
  pot.omeg.loc <- rescaled.posterior.potentials[,3]

  pot.tau1.int <- HPDI(samples = pot.tau1.loc, prob = prob.level)
  pot.tau2.int <- HPDI(samples = pot.tau2.loc, prob = prob.level)
  pot.omeg.int <- HPDI(samples = pot.omeg.loc, prob = prob.level)
  pot.ints     <- rbind(pot.tau1.int, pot.tau2.int, pot.omeg.int)

  pot.indics     <- !cbind( pot.ints[,1] <= 1 & pot.ints[,2] >= 1 )
  pot.indics.txt <- c("","","")
  if(pot.indics[1] == TRUE){
    pot.indics.txt[1] <- paste0(loc.node.names[1]," not 50/50")    # Reject X1 1 and X1 2 about the same frequency
  } else {
    pot.indics.txt[1] <- paste0(loc.node.names[1]," around 50/50") # Don't reject X1 1 and X1 2 about the same frequency
  }
  if(pot.indics[2] == TRUE){
    pot.indics.txt[2] <- paste0(loc.node.names[2]," not 50/50")    # Reject X2 1 and X2 2 about the same frequency
  } else {
    pot.indics.txt[2] <- paste0(loc.node.names[2]," around 50/50") # Don't reject X2 1 and X2 2 about the same frequency
  }
  if(pot.indics[3] == TRUE){
    pot.indics.txt[3] <- paste0(loc.node.names[1]," NOT independent ",loc.node.names[2]) # Reject X1 _||_ X2
  } else {
    pot.indics.txt[3] <- paste0(loc.node.names[1]," independent ",loc.node.names[2])     # Don't reject X1 _||_ X2
  }

  pot.ints     <- data.frame( pot.ints, pot.indics, pot.indics.txt )
  rownames(pot.ints) <- c("e^tau1", "e^tau2", "e^omega")
  colnames(pot.ints) <- c(paste0("|",prob.level), paste0(prob.level,"|"), "1.NOT.covered?", "Indication")

  if(printQ==TRUE){
    print("=========== Indicators with respect to potentials ===========")
    print(pot.ints)
    print("=============================================================")
  }

  # Posterior energies
  tau1.loc     <- log(rescaled.posterior.potentials[,1])
  tau2.loc     <- log(rescaled.posterior.potentials[,2])
  omeg.loc     <- log(rescaled.posterior.potentials[,3])

  tau1.int   <- HPDI(samples = tau1.loc, prob = prob.level)
  tau2.int   <- HPDI(samples = tau2.loc, prob = prob.level)
  omeg.int   <- HPDI(samples = omeg.loc, prob = prob.level)
  ene.ints     <- rbind(tau1.int, tau2.int, omeg.int)

  ene.indics     <- !cbind( ene.ints[,1] <= 0 & ene.ints[,2] >= 0 )
  ene.indics.txt <- c("","","")
  if(ene.indics[1] == TRUE){
    ene.indics.txt[1] <- paste0(loc.node.names[1]," not 50/50")    # Reject X1 1 and X1 2 about the same frequency
  } else {
    ene.indics.txt[1] <- paste0(loc.node.names[1]," around 50/50") # Don't reject X1 1 and X1 2 about the same frequency
  }
  if(ene.indics[2] == TRUE){
    ene.indics.txt[2] <- paste0(loc.node.names[2]," not 50/50")    # Reject X2 1 and X2 2 about the same frequency
  } else {
    ene.indics.txt[2] <- paste0(loc.node.names[2]," around 50/50") # Don't reject X2 1 and X2 2 about the same frequency
  }
  if(ene.indics[3] == TRUE){
    ene.indics.txt[3] <- paste0(loc.node.names[1]," NOT independent ",loc.node.names[2]) # Reject X1 independent X2
  } else {
    ene.indics.txt[3] <- paste0(loc.node.names[1]," independent ",loc.node.names[2])     # Don't reject X1 independent X2
  }

  ene.ints     <- data.frame(ene.ints, ene.indics, ene.indics.txt )
  rownames(ene.ints) <- c("tau1", "tau2", "omega")
  colnames(ene.ints) <- c(paste0("|",prob.level), paste0(prob.level,"|"), "0.NOT.covered?", "Indication")

  if(printQ==TRUE) {
    print("============ Indicators with respect to energies ============")
    print(ene.ints)
    print("=============================================================")
  }

  marginal.edge.interval.indicator.info <- list(
    pot.ints,
    ene.ints
  )
  names(marginal.edge.interval.indicator.info) <- c("potentials.int.info","energies.int.info")

  return(marginal.edge.interval.indicator.info)

}


#' marginal edge bayes independence indicators from beliefs
#'
#' Outputs statistics on posteriors for Pr(X1), Pr(X1|X2), Pr(X2), Pr(X2|X1)
#' Handy for checking for independencies between X1 and X2
#'
#' @param XX XX
#'
#' @return The function will XX
#'
#'
#' @export
marginal.edge.bayes.beleifs.indep.indic <- function(posterior.edge.belief.info, prob.level=0.95, printQ=FALSE){

  distc.nmes <- names(posterior.edge.belief.info)[c(4,5)] # Edge conditional distributions
  distm.nmes <- names(posterior.edge.belief.info)[c(2,3)] # Edge marginal node distributions

  sect.txts       <- array("", c(8,2))
  sect.indic.vals <- array(logical(0), c(8,1))
  sect.cond.pis   <- array(NA, c(8,2))
  sect.marg.pis   <- array(NA, c(8,2))
  sect.marg.meds  <- array(NA, c(8,1))

  # i1 = X1|X2, i2 = X2|X1
  count <- 1
  for(i in 1:2) {

    distc.nme      <- distc.nmes[i]
    distm.nme      <- distm.nmes[i]
    cond.sect.nmes <- colnames(posterior.edge.belief.info[[distc.nme]])
    marg.sect.nmes <- colnames(posterior.edge.belief.info[[distm.nme]])

    # sector 1,1 -> Pr(XA=1 | XB=1), Pr(XA = 1)
    cond.sect.nme <- cond.sect.nmes[1]
    marg.sect.nme <- marg.sect.nmes[1]
    med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
    cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
    marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)
    med.in.intQ   <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
    sect.cond.pis[count, ] <- cond.pi
    sect.marg.meds[count]  <- med.lne
    sect.marg.pis[count, ] <- marg.pi
    sect.txt      <- paste0("Bel(",cond.sect.nme,") == Bel(",marg.sect.nme,")")
    if(med.in.intQ == TRUE) {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "YES"
      sect.indic.vals[count] <- TRUE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    } else {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "NO"
      sect.indic.vals[count] <- FALSE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    }
    count <- count + 1

    # sector 2,1 -> Pr(XA=2 | XB=1), Pr(XA = 2)
    cond.sect.nme <- cond.sect.nmes[3]
    marg.sect.nme <- marg.sect.nmes[2]
    med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
    cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
    marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)
    med.in.intQ   <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
    sect.cond.pis[count, ] <- cond.pi
    sect.marg.meds[count]  <- med.lne
    sect.marg.pis[count, ] <- marg.pi
    sect.txt      <- paste0("Bel(",cond.sect.nme,") == Bel(",marg.sect.nme,")")
    if(med.in.intQ == TRUE) {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "YES"
      sect.indic.vals[count] <- TRUE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    } else {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "NO"
      sect.indic.vals[count] <- FALSE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    }
    count <- count + 1

    # sector 1,2 -> Pr(XA=1 | XB=2), Pr(XA = 1)
    cond.sect.nme <- cond.sect.nmes[2]
    marg.sect.nme <- marg.sect.nmes[1]
    med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
    cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
    marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)
    med.in.intQ   <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
    sect.cond.pis[count, ] <- cond.pi
    sect.marg.meds[count]  <- med.lne
    sect.marg.pis[count, ] <- marg.pi
    sect.txt      <- paste0("Bel(",cond.sect.nme,") == Bel(",marg.sect.nme,")")
    if(med.in.intQ == TRUE) {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "YES"
      sect.indic.vals[count] <- TRUE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    } else {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "NO"
      sect.indic.vals[count] <- FALSE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    }
    count <- count + 1

    # sector 2,2 -> Pr(XA=2 | XB=2), Pr(XA = 2)
    cond.sect.nme <- cond.sect.nmes[4]
    marg.sect.nme <- marg.sect.nmes[2]
    med.lne       <- median(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme])
    cond.pi       <- HPDI(posterior.edge.belief.info[[distc.nme]][,cond.sect.nme], prob.level)
    marg.pi       <- HPDI(posterior.edge.belief.info[[distm.nme]][,marg.sect.nme], prob.level)
    med.in.intQ   <- (cond.pi[1] <= med.lne & cond.pi[2] >= med.lne)
    sect.cond.pis[count, ] <- cond.pi
    sect.marg.meds[count]  <- med.lne
    sect.marg.pis[count, ] <- marg.pi
    sect.txt      <- paste0("Bel(",cond.sect.nme,") == Bel(",marg.sect.nme,")")
    if(med.in.intQ == TRUE) {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "YES"
      sect.indic.vals[count] <- TRUE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    } else {
      sect.txts[count,1]     <- sect.txt
      sect.txts[count,2]     <- "NO"
      sect.indic.vals[count] <- FALSE
      #print(paste0(sect.txt, " ==> ", sect.txts[count,2]))
    }
    count <- count + 1

  }

  meidp.info <- data.frame(sect.txts, sect.indic.vals, sect.cond.pis, sect.marg.meds, sect.marg.pis)
  colnames(meidp.info) <- c(
    "H0",
    "Accept H0?",
    "H0.Q",
    paste0("cond: |",prob.level), paste0(prob.level, "|"),
    "Med[Bel(XA)]",
    paste0("marg: |",prob.level), paste0(prob.level, "|")
  )

  return(meidp.info)

}
