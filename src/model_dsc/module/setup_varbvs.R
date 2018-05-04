X <- data$X
storage.mode(X) <- "double"
n <- nrow(X)
p <- ncol(X)
X <- scale(X,center = TRUE,scale = FALSE)
alpha0  <- runif(p)
alpha0  <- alpha0/sum(alpha0)
mu0     <- rnorm(p)
pp      <- rep(maxL/p, p)
logodds <- varbvs:::logit(pp)
Y <- data$Y
for (r in 1:ncol(Y)) {
  Y[,r] <- Y[,r] - mean(Y[,r])
}
storage.mode(Y) <- "double"
