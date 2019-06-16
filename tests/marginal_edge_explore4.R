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

infer.exact(junk)

