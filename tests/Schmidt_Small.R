library(CRFutil)

#The potentials for Cathy-Heather-Mark-Allison: 1-2-3-4
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi23 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi34 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))

# Define states and feature function:
s1 <- "right"
s2 <- "wrong"
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("Cathy", "Heather", "Mark", "Alison")

edges <- rbind(
  c(1,2), #Cathy-Heather
  c(2,3), #Heather-Mark
  c(3,4)  #Mark-Allison
)

# Define a convenience function wrapper:
ener.func <- function(config) {
  engy <-config.energy(config = config, edges.mat = edges,
                       one.lgp = list(log(Psi1),log(Psi2),log(Psi3),log(Psi4)),
                       two.lgp = list(log(Psi12),log(Psi23),log(Psi34)), ff = f)
  return(engy)
}

# All configuration energies:
config.energies <- sapply(1:nrow(config.mat), function(xx){ener.func(config.mat[xx,])})
prodPots        <- exp(config.energies)
Z               <- sum(prodPots)
Prs             <-prodPots/Z

Z
sum(Prs)
cbind(config.mat, config.energies, prodPots, Prs)

