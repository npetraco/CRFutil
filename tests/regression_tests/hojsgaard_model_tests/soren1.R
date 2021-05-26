library(gRbase)
library(gRim)

data(reinis)
class(reinis)
dim(reinis)   # a 2^6 table
reinis

grphf <- ~smoke:systol + smoke:mental:phys
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

dm1 <- dmod(formula = grphf, data = reinis, interactions = 3, fit = T)
dm1
summary(dm1)
terms(dm1)
formula(dm1)
dm1$glist
class(dm1)
print(dm1)
logLik(dm1)
modelProperties(dm1)

?loglin
dmod
dm1$fitinfo$dimension

reinisAGG <- as.data.frame(reinis)
reinisAGG

attach(reinisAGG)
gm1 <- glm(Freq~-1+smoke*systol + smoke*mental*phys,family=poisson, data=reinisAGG)
#gm1 <- glm(Freq~-1+smoke:systol + smoke:mental:phys,family=poisson, data=reinisAGG)
summary(gm1)
detach(reinisAGG)

# smoke mental phys systol protein family
im1 <- loglin(table = reinis, margin = list(c(1,4),c(1,2,3)), param=T)
im1$pearson
dm1$fitinfo$pearson

im1$param

library(MASS)
?loglm
minn38
xtabs(Freq~., reinisAGG)
im2 <- loglm(~smoke*systol + smoke*mental*phys - 1, reinis)
im2$param
