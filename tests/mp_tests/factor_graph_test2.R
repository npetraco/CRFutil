library(CRFutil)

# Model:
grphf <- ~1:2 + 1:3 + 1:4 + 1:5

# Convert to pw factor graph
pwfg <- mrf2pwfg(grphf, plotQ=T)

# Schedules:
schs <- get.root.paths(pwfg, root.node = 5, serial.schedsQ = T)
schs$forwrd

# Make up some potentials for testing:
adj <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate a model:
knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 9)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$node.par[4,1,] <- 4
knm$node.par[5,1,] <- 5
knm$edge.par[[1]][1,1,1] <- 6
knm$edge.par[[1]][2,2,1] <- 6
knm$edge.par[[2]][1,1,1] <- 7
knm$edge.par[[2]][2,2,1] <- 7
knm$edge.par[[3]][1,1,1] <- 8
knm$edge.par[[3]][2,2,1] <- 8
knm$edge.par[[4]][1,1,1] <- 9
knm$edge.par[[4]][2,2,1] <- 9

#set.seed(1)
knm$par <- runif(9,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
#num.samps <- 1000
#set.seed(1)
#samps <- sample.exact(knm, num.samps)
#mrf.sample.plot(samps)

out.pot
gr.pots <- make.gRbase.potentials(crf = knm, node.names = as.character(1:5), state.nmes = c("up", "dn"))
ep <- gr.pots$edge.potentials
np <- gr.pots$node.potentials
np
ep

ep[[1]] %a*% ep[[2]]
junk <- ep[[1]] %a*% ar_prod_list(list(ep[[1]], ep[[2]], ep[[3]]))
dim(junk)
names(dimnames(junk))
ar_marg(ep[[1]],2)

schs$forwrd
make.f2v.msg(
  in.v.msgs.list= lapply(c(1,3), function(xx){ep[[xx]]}),
  f.msg = ep[[4]],
  out.v.nme = 1)

schs$forwrd
m.f1.1 <- make.f2v.msg(f.msg = np[[1]], out.v.nme = 1)
m.f2.2 <- make.f2v.msg(f.msg = np[[2]], out.v.nme = 2)
m.f3.3 <- make.f2v.msg(f.msg = np[[3]], out.v.nme = 3)
m.f4.4 <- make.f2v.msg(f.msg = np[[4]], out.v.nme = 4)
m.f5.5 <- make.f2v.msg(f.msg = np[[5]], out.v.nme = 5)

m.1.f15 <- make.v2f.msg(list(m.f1.1))
m.2.f12 <- make.v2f.msg(list(m.f2.2))
m.3.f13 <- make.v2f.msg(list(m.f3.3))
m.4.f14 <- make.v2f.msg(list(m.f4.4))

m.f15.5 <- make.f2v.msg(in.v.msgs.list = list(m.1.f15), f.msg = ep[[4]], out.v.nme = 5)


