library(CRFutil)

# Slayer field model:
grphf <- ~1:2+1:3+1:4+1:5+2:3+2:4+2:5+3:4+3:5
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)


# Make up random potentials and return a CRF-object
slay    <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=25, seed=1)
samps       <- slay$samples
known.model <- slay$model
known.model$node.par
    known.model$edges[9,]
known.model$edge.par[[9]]


# Feature function
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

# Grab the first sampled configuration:
X <- samps[1,]
X

known.model$par

nlpl.info <- neglogpseudolik.config(
  param        = NULL,
  config       = X,
  crf          = known.model,
  cond.en.form = "feature",
  gradQ        = T,
  ff           = f0)
nlpl.info


gnlpl.info <-grad.neglogpseudolik.config(
  param                      = NULL,
  config                     = X,
  crf                        = known.model,
  cond.en.form               = "feature",
  ff                         = f0)

gnlpl.info$grad.neglogpseudolik


samp.lik.info <- neglogpseudolik(
  par         = NULL,
  crf         = known.model,
  samples     = samps,
  conditional.energy.function.type="feature",
  ff          = f0,
  gradQ       = T,
  update.crfQ = F)

samp.lik.info$samp.neglogpseudolik
sum(samp.lik.info$nliks)
samp.lik.info$samp.grad.neglogpseudolik
colSums(samp.lik.info$gnliks)

