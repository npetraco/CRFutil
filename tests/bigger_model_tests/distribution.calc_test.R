m1 <- distribution.calc(fit)
m1$logZ

m2 <- distribution.calc(fit, logZ.calc = infer.exact)
m2$logZ

m3 <- distribution.calc(fit, logZ.calc = infer.junction)
m3$logZ

m4 <- distribution.calc(fit, logZ.calc = infer.lbp)
m4$logZ

m5 <- distribution.calc(fit, logZ.calc = infer.rbp)
m5$logZ

m6 <- distribution.calc(fit, logZ.calc = infer.trbp)
m6$logZ

m7 <- distribution.calc(fit, logZ.calc = infer.tree)
m7$logZ

m8 <- distribution.calc(fit, logZ.calc = infer.cutset)
m8$logZ


logZ.mle.model

KLD(joint.mle[,pr.idx], m1$state.probs)
KLD(joint.mle[,pr.idx], m2$state.probs)
KLD(joint.mle[,pr.idx], m3$state.probs)
KLD(joint.mle[,pr.idx], m4$state.probs)
KLD(joint.mle[,pr.idx], m5$state.probs)
KLD(joint.mle[,pr.idx], m6$state.probs)
KLD(joint.mle[,pr.idx], m7$state.probs)
KLD(joint.mle[,pr.idx], m8$state.probs)



class(infer.exact)
is.null(infer.exact)
