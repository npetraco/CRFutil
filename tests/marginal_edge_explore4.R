library(CRFutil)

# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 100000

# Must label sample states 1 and 2
samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.5)
)


mrf.sample.plot(samps)
#samps
junk.e <- marginal.edge.emp.pr(samps, printQ = T)

junk <- NULL
junk <- marginal.edge.mrf(samps)
dump.crf(crf = junk)

junk.be <- infer.exact(junk)
junk.be$node.bel
junk.b1

junk$edges
make.gRbase.beliefs(junk.be, node.names = c("X1", "X2"), junk$edges, state.nmes=c("1","2"))

marginal.edge.bels(edge.mrf.obj = junk, node.names = c("A","B"), printQ = T)
