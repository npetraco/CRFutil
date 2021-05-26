library(CRFutil)
library(gRbase)
library(Rgraphviz)
library(MASS)
library(gRim)

grf.eq <- ~A + B + A:B
AB <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)
AB$edge.pot[[1]][1,1] <- AB$edge.pot[[1]][2,2] <- 10 # ******* theta true 11,22 is theta-times more likely than 12 or 21

num.samps <- 100
samps     <- sample.exact(AB, num.samps)
mrf.sample.plot(samps)
#samps

c("diam","height","species")
colnames(samps) <- c("diam","height")
head(samps)

X  <- samps                              # Raw Observed states
Xc <- xtabs(~., data=data.frame(X))                  # Contingency table
Xf <- as.data.frame(xtabs(~., data=data.frame(X)))   # Freq table of observed states
X
Xc
Xf

# Model matrix ??????
Xf
Xf[,1]
dia.dc <- Xf[,1]
contrasts(dia.dc)
contrasts(dia.dc) <- contr.sum(2)
contrasts(dia.dc)

Xf[,2]
hgt.dc <- Xf[,2]
contrasts(hgt.dc)
contrasts(hgt.dc) <- contr.sum(2)
contrasts(hgt.dc)

# Xf[,3]
# spc.dc <- Xf[,3]
# contrasts(spc.dc)
# contrasts(spc.dc) <- contr.sum(2)
# contrasts(spc.dc)

Xm <- model.matrix(~dia.dc + hgt.dc + dia.dc:hgt.dc,
                   contrasts = list(dia.dc = "contr.sum",
                                    hgt.dc = "contr.sum"))
Xm

loglm( ~ diam + height + diam:height, data = Xc)
loglm(Freq~dia.dc + hgt.dc + dia.dc:hgt.dc, data = Xf)

c2 <- coef(loglm( ~ diam + height + diam:height, data = Xc))
c3 <- coef(loglm(Freq~dia.dc + hgt.dc + dia.dc:hgt.dc, data = Xf))

c2[[1]]
c3[[1]]

c2[[2]]
c3[[2]] # ??

c2[[3]]
c3[[3]] # ??

c2[[4]]
c3[[4]] # ??

# glm next
glm.f  <- glm(Freq~dia.dc + hgt.dc + dia.dc:hgt.dc, data = Xf, family = poisson(link="log"))
#glm.f2 <- glm(Freq~diam + height + diam:height, data = Xf, family = poisson(link="log"))

lmcoefs <- c(c2[[1]][1], c2[[2]][1], c2[[3]][1], c2[[4]][1]) # loglm coefs
lmcoefs
names(lmcoefs) <- names(c2)

lmcoefs     # compare
coef(glm.f)
coef(glm.f2)

summary(glm.f)
#summary(glm.f2)

# According to p-values pf coefs:
par.glm <- c(0,0,0.86560)

grf.eq <- ~A+B+A:B
AB <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)
AB$node.pot[1,1] <- exp(par.glm[1])
AB$node.pot[2,1] <- exp(par.glm[2])
AB$node.pot[1,2] <- exp(-par.glm[1])
AB$node.pot[2,2] <- exp(-par.glm[2])

AB$edge.pot[[1]][1,1] <- AB$edge.pot[[1]][2,2] <- exp(par.glm[3])
AB$edge.pot[[1]][1,2] <- AB$edge.pot[[1]][2,1] <- exp(-par.glm[3])

AB$node.pot
AB$edge.pot

# Re-scale?? Get about theta true??
AB$edge.pot[[1]][1,]/min(AB$edge.pot[[1]][1,])
AB$edge.pot[[1]][2,]/min(AB$edge.pot[[1]][2,])

marginal.edge.loglin(samps)
samps
xtabs(~., data=data.frame(samps))
