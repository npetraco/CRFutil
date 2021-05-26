x <- seq(from = -2.0, to=2.0, by = 0.05)
y <- seq(from = -2.0, to=2.0, by = 0.05)
xy.grid <- as.matrix(expand.grid(x,y))
dim(xy.grid)

# function:
f <- function(xx,yy){
  xx*exp(-(xx^2 + yy^2))
}
# gradient:
gf <- function(xx,yy){
  c(exp(-xx^2 - yy^2) - 2*exp(-xx^2 - yy^2)*xx^2,-2*exp(-xx^2 - yy^2)*xx*yy)
}

# Gradient of the sample neg. log pseudo-likelihood:
# Compute and plot gradient:
f.mat        <- array(NA, c(length(x), length(y)))
grad.f.mat.u <- array(NA, c(length(x), length(y)))
grad.f.mat.v <- array(NA, c(length(x), length(y)))

count <- 1
for(i in 1:length(x)) {
  for(j in 1:length(y)) {


    print(count)
    f.mat[i,j] <- f(x[i], y[j])
    gfij       <- gf(x[i], y[j])

    grad.f.mat.u[i,j] <- gfij[1]
    grad.f.mat.v[i,j] <- gfij[2]

    count <- count + 1
  }
}

hist(f.mat)
hist(grad.f.mat.u)
hist(grad.f.mat.v)

# Make a 2D plot of the gradient ???????:
library(OceanView)
u <- grad.f.mat.u
v <- grad.f.mat.v

# Sample negative log pseudo-likelihood:
image2D(x = x, y = y, z = f.mat, contour = TRUE, xlab="x", ylab="y")

# Gradient:
quiver2D(x = x, y = y, u = u, v = v, add = TRUE, by = 6)

