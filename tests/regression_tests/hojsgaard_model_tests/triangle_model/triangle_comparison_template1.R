# Compare the Bayes, and MLE with the True model:
bayes.lrm$par <- c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,4]), median(beta[,5]), median(beta[,6]))
mle.lrm$par   <- as.numeric(coef(Ma))

out.pot2  <- make.pots(parms = bayes.lrm$par, crf = bayes.lrm, rescaleQ = T, replaceQ = T)
out.pot2a <- make.pots(parms = mle.lrm$par,   crf = mle.lrm,   rescaleQ = T, replaceQ = T)

# Node and edge beliefs:
bayes.lrm.bel  <- infer.exact(bayes.lrm)
mle.lrm.bel    <- infer.exact(mle.lrm)
true.model.bel <- infer.exact(true.model)

bayes.lrm.bel$node.bel
mle.lrm.bel$node.bel
true.model.bel$node.bel

bayes.lrm.bel$edge.bel[[1]]
mle.lrm.bel$edge.bel[[1]]
true.model.bel$edge.bel[[1]]

bayes.lrm.bel$edge.bel[[2]]
mle.lrm.bel$edge.bel[[2]]
true.model.bel$edge.bel[[2]]

bayes.lrm.bel$edge.bel[[3]]
mle.lrm.bel$edge.bel[[3]]
true.model.bel$edge.bel[[3]]

# True configuration probabilities:
pot.info.true.model       <- make.gRbase.potentials(true.model, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials, pot.info.true.model$edge.potentials)
logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))

# Bayes configuration probabilities:
pot.info        <- make.gRbase.potentials(bayes.lrm, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))

# MLE configuration probabilities:
pot.info2        <- make.gRbase.potentials(mle.lrm, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info2    <- distribution.from.potentials(pot.info2$node.potentials, pot.info2$edge.potentials)
logZ2            <- gR.dist.info2$logZ
joint.dist.info2 <- as.data.frame(as.table(gR.dist.info2$state.probs))

# Compare logistic regression based config probs from Bayes and MLE
# with the true config probs:
bayes.lrm.cp   <- round(joint.dist.info[,4]*100, 1)            # Bayes logistic
mle.lrm.cp     <- round(joint.dist.info2[,4]*100, 1)           # MLE logistic
true.model.cp  <- round(joint.dist.info.true.model[,4]*100, 1) # True
cbind(joint.dist.info[,c(2,3,1)], bayes.lrm.cp, mle.lrm.cp, true.model.cp)
