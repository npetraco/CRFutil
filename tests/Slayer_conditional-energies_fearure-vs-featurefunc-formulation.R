library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")

# Make up random potentials and return a CRF-object
num.samps   <- 100
n.states    <- 2
slay        <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states,
                                num.sims=num.samps, seed=1)
samps       <- slay$samples
known.model <- slay$model
mrf.sample.plot(samps)

# Needed for the energy functions:
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function


# For storage:
en.result <- array(0, c(num.samps*ncol(samps),5))

# Loop around the sample numbers and then the elements of each sample
count <- 1
for(i in 1:num.samps) {
  for(j in 1:ncol(samps)) {
    samp.num <- i
    elem.num <- j
    # Feature formulation:
    ce1 <- conditional.config.energy(
      config = samps[samp.num,],
      condition.element.number = elem.num,
      crf=known.model,
      ff=f0, printQ=FALSE)

    # Feature function formulation:
    ce2 <- conditional.config.energy2(
             config = samps[samp.num,],
             condition.element.number = elem.num,
             crf=known.model,
             ff=f0, printQ=FALSE)

    en.result[count,] <- c(i,j,ce1,ce2,ce1-ce2)
    count <- count + 1
  }
}
en.result
en.result[,5]
plot(1:length(en.result[,5]),en.result[,5], ylab="Difference", xlab="Config index")
