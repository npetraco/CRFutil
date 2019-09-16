
# Load results
setwd("/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/tesseract_model/")

load(file = "tesseract_mle_model.dist.RData")
mle.info <- joint.dist.info.ordered

load(file = "tesseract_bayes_loglin_model.dist.RData")
llm.info <- joint.dist.info.ordered

load(file = "tesseract_true_model.dist.RData")
tru.info <- joint.dist.info.ordered

config.prs.mle <- mle.info[, ncol(mle.info)]
config.prs.llm <- llm.info[, ncol(llm.info)]
config.prs.tru <- tru.info[, ncol(tru.info)]


plot(config.prs.mle*100, typ="h", ylim = c(0,0.65))
plot(config.prs.llm*100, typ="h", ylim = c(0,0.65))
plot(config.prs.tru*100, typ="h", ylim = c(0,0.65))

tru.v.mle <- config.prs.tru - config.prs.mle
tru.v.llm <- config.prs.tru - config.prs.llm
plot(tru.v.mle*100, typ="h", main="Fit diffs (%)")
plot(tru.v.llm*100, typ="h", main="Fit diffs (%)")

hist(tru.v.mle*100)
hist(tru.v.llm*100)

par( mfrow = c( 2, 1 ) )
plot(config.prs.mle*100, typ="h", ylim = c(0,0.65))
plot(tru.v.mle*100, typ="h", main="Fit diffs (%)")
dev.off()

par( mfrow = c( 2, 1 ) )
plot(config.prs.llm*100, typ="h", ylim = c(0,0.65))
plot(tru.v.llm*100, typ="h", main="Fit diffs (%)")
dev.off()

par( mfrow = c( 2, 1 ) )
plot(config.prs.mle*100, typ="h", ylim = c(0,0.65))
plot(config.prs.tru*100, typ="h", ylim = c(0,0.65))
dev.off()

par( mfrow = c( 2, 1 ) )
plot(config.prs.llm*100, typ="h", ylim = c(0,0.65))
plot(config.prs.tru*100, typ="h", ylim = c(0,0.65))
dev.off()
