library(magrittr)
library(ISwR)
library(rstan)

# Danish smoker data
data(eba1977)
summary(eba1977)
levels(eba1977[,1])


glm1 <- glm(formula = cases ~ age + city + offset(log(pop)),
            family  = poisson(link = "log"),
            data    = eba1977)
summary(glm1)

## Model matrix
modMat <- as.data.frame(model.matrix(glm1))
modMat$offset <- log(eba1977$pop)
names(modMat) <- c("intercept", "age55_59", "age60_64", "age65_69", "age70_74",
                   "age75plus", "cityHorsens", "cityKolding", "cityVejle", "offset")
modMat

dat   <- as.list(modMat)
dat$y <- eba1977$cases
dat$N <- nrow(modMat)
dat$p <- ncol(modMat) - 1

## Load Stan file
fileName <- "tests/regression_tests/hojsgaard_model_tests/stan_loglin_test1.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
cat(stan_code)

## Run Stan
resStan <- stan(model_code = stan_code, data = dat, chains = 3, iter = 3000, warmup = 500, thin = 10)
