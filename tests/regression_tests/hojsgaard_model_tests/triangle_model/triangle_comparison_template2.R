
# loglin gRim,               DONE
# bayes (poisson regr) loglin
# zero inflated loglin
# bayes zero inflated loglin
# model matrix loglin,       DONE
# glm-poisson,               DONE, see triangle-model_loglin-modelmartix, bottom. Equivalent results to model matrix loglin and loglin gRim
# bayes logistic,            DONE
# mle logistic,              DONE
# CRF fit                    DONE
# PSL fit


library(CRFutil)
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/triangle_data/"
#setwd(fpth)

load(paste0(fpth,"triangle_loglin_dist.RData"))
loglin.dist.info

load(paste0(fpth,"triangle_exact_dist.RData"))
exact.dist.info

load(paste0(fpth,"triangle_empirical_dist.RData"))
empirical.dist.info

load(paste0(fpth,"triangle_loglin_MM_dist.RData"))
loglin.MM.dist.info

load(paste0(fpth,"triangle_bayes-lr_dist.RData"))
bayes.lr.dist.info
bayes.lr.dist.info <- bayes.lr.dist.info[,c(2,3,1,4)] # Forgot to rearrange columns

load(paste0(fpth,"triangle_mle-lr_dist.RData"))
mle.lr.dist.info

load(paste0(fpth,"triangle_mrf_dist.RData"))
mrf.dist.info


ref.states <- exact.dist.info[,1:3]
loglin.rearr.idxs    <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = loglin.dist.info[,1:3])})
emp.rearr.idxs       <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = empirical.dist.info[,1:3])})
loglin.MM.rearr.idxs <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = loglin.MM.dist.info[,1:3])})
bayes.lr.rearr.idxs  <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = bayes.lr.dist.info[,1:3])})
mle.lr.rearr.idxs    <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = mle.lr.dist.info[,1:3])})
mrf.rearr.idxs       <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = mrf.dist.info[,1:3])})


loglin.rearr.idxs
loglin.MM.rearr.idxs
emp.rearr.idxs
bayes.lr.rearr.idxs
mle.lr.rearr.idxs
mrf.rearr.idxs


spcr <- rep(NA, nrow(ref.states))
cbind(exact.dist.info[,1:3], spcr,
      loglin.dist.info[loglin.rearr.idxs,1:3], spcr,
      loglin.MM.dist.info[loglin.MM.rearr.idxs,1:3], spcr,
      bayes.lr.dist.info[bayes.lr.rearr.idxs,1:3], spcr,
      mle.lr.dist.info[mle.lr.rearr.idxs,1:3], spcr,
      mrf.dist.info[mle.lr.rearr.idxs,1:3], spcr,
      empirical.dist.info[emp.rearr.idxs,1:3])

loglin.dist.info    <- loglin.dist.info[loglin.rearr.idxs,]
loglin.MM.dist.info <- loglin.MM.dist.info[loglin.MM.rearr.idxs,]
empirical.dist.info <- empirical.dist.info[emp.rearr.idxs,]
bayes.lr.dist.info  <- bayes.lr.dist.info[bayes.lr.rearr.idxs,]
mle.lr.dist.info    <- mle.lr.dist.info[mle.lr.rearr.idxs,]
mrf.dist.info       <- mrf.dist.info[mrf.rearr.idxs,]


cbind(exact.dist.info[,1:3], spcr,
      loglin.dist.info[,1:3], spcr,
      loglin.MM.dist.info[,1:3], spcr,
      bayes.lr.dist.info[,1:3], spcr,
      mle.lr.dist.info[,1:3], spcr,
      mrf.dist.info[,1:3], spcr,
      empirical.dist.info[,1:3])

cbind(ref.states,
      exact.dist.info[,4],
      loglin.dist.info[,4],
      loglin.MM.dist.info[,4],
      bayes.lr.dist.info[,4],
      mle.lr.dist.info[,4],
      mrf.dist.info[,4],
      empirical.dist.info[,4])
