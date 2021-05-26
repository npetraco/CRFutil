library(CRFutil)

# Model:
grphf <- ~1:2 # Node names must be numbers for now.... FIX

# Convert to pw factor graph
pwfg <- mrf2pwfg(grphf, plotQ=T)

# Schedules:
schs <- get.root.paths(pwfg, root.node = 2, serial.schedsQ = T)
schs$forwrd
schs$backward

# Make up some potentials for testing:
adj <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate a model:
knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 3)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$edge.par[[1]][1,1,1] <- 3
knm$edge.par[[1]][2,2,1] <- 3

#set.seed(1)
knm$par <- runif(3,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
#num.samps <- 1000
#set.seed(1)
#samps <- sample.exact(knm, num.samps)
#mrf.sample.plot(samps)

out.pot
gr.pots <- make.gRbase.potentials(crf = knm, node.names = as.character(1:5), state.nmes = c("up", "dn"))
ep <- gr.pots$edge.potentials # Edge potentials in gRbase format
np <- gr.pots$node.potentials # Node potentials in gRbase format
np
ep

# Normalizes node potentials
class(np)
class(np[[1]])
sapply(1:length(np), function(xx){tabNormalize(tab = np[[xx]], type = 1)})

# Normalize edge potentials
tabNormalize(tab = ep[[1]], type = 1) # Normalizes rows of edge potential
tabNormalize(tab = ep[[1]], type = 2) # Normalizes whole edge potential

# Marginalization
tabMarg(tab = ep[[1]], marg = c(1)) # Marginalize out all variables except 1
tabMarg(tab = ep[[1]], marg = c(2)) # Marginalize out all variables except 2
tabMarg(tab = ep[[1]])              # Sum all edge potential entries
tabMarg(tab = np[[1]])              # Sum all node 1 potential entries
tabMarg(tab = np[[2]])              # Sum all node 2 potential entries

tabPerm(ep[[1]], perm = c(2,1))     # Permute order of variable defining the edge potential

np
tabMult(np[[1]], np[[2]])
tabMult(ep[[1]], np[[1]])
tabMult(ep[[1]], np[[2]])
tabMult(ep[[1]], ep[[1]])

# Need neighbor function with ability to drop some neighbors, or we have already in makef,v??
igpwfg <- graph_from_graphnel(pwfg)
V(igpwfg)
neighborhood(igpwfg, nodes = V(igpwfg)[4], mindist = 1)





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


