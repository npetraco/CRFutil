
# loglin gRim,         DONE
# model matrix loglin, DONE
# glm-poisson,         DONE, see triangle-model_loglin-modelmartix, bottom. Equivalent results to model matrix loglin
# bayes logistic,
# mle logistic,
# CRF fit
# PSL fit


library(CRFutil)
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/triangle_data/"

load(paste0(fpth,"triangle_loglin_dist.RData"))
loglin.dist.info

load(paste0(fpth,"triangle_exact_dist.RData"))
exact.dist.info

load(paste0(fpth,"triangle_empirical_dist.RData"))
empirical.dist.info

load(paste0(fpth,"triangle_loglin_MM_dist.RData"))
loglin.MM.dist.info


ref.states <- exact.dist.info[,1:3]
loglin.rearr.idxs    <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = loglin.dist.info[,1:3])})
emp.rearr.idxs       <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = empirical.dist.info[,1:3])})
loglin.MM.rearr.idxs <- sapply(1:nrow(ref.states), function(xx){row.match(ref.states[xx,], table = loglin.MM.dist.info[,1:3])})


loglin.rearr.idxs
loglin.MM.rearr.idxs
emp.rearr.idxs

spcr <- rep(NA, nrow(ref.states))
cbind(exact.dist.info[,1:3], spcr,
      loglin.dist.info[loglin.rearr.idxs,1:3], spcr,
      loglin.MM.dist.info[loglin.MM.rearr.idxs,1:3], spcr,
      empirical.dist.info[emp.rearr.idxs,1:3])


loglin.dist.info    <- loglin.dist.info[loglin.rearr.idxs,]
loglin.MM.dist.info <- loglin.MM.dist.info[loglin.MM.rearr.idxs,]
empirical.dist.info <- empirical.dist.info[emp.rearr.idxs,]


cbind(exact.dist.info[,1:3], spcr,
      loglin.dist.info[,1:3], spcr,
      loglin.MM.dist.info[,1:3], spcr,
      empirical.dist.info[,1:3])

cbind(ref.states,
      exact.dist.info[,4],
      loglin.dist.info[,4],
      loglin.MM.dist.info[,4],
      empirical.dist.info[,4])
