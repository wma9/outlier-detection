#####################
### simulate data ###
##### Settings ######
#####################

N <- 100
width <- 25
b <- 5
c <- 5
sigma2_0 <- 2
sigma2_c <- 20
alpha <- 0.05
K <- 5 #question

library(fields)
library(geoR)
set.seed(123)
lat <- runif(N, 0, 1)*width
lon <- runif(N, 0, 1)*width
s <- cbind(lat, lon)
beta <- rnorm(N, 0, 1)
d <- as.matrix(dist(s))
covariance <- cov.spatial(d, cov.model= "exponential", cov.pars=c(b, c))
# covariance <- apply(covariance, 2, function(x){ifelse(x<b*exp(-1), 0, x)})
w <- t(chol(covariance))%*%rnorm(N)
eps <- rep(0, N)
for (i in 1:N){
  test <- runif(1)
  eps[i] <- ifelse(test > alpha, rnorm(1, 0, sqrt(sigma2_0)), rnorm(1, 0, sqrt(sigma2_c)))
}

Z <- beta + w + eps

Y <- cbind(Z, lat, lon)
Y <- data.frame(dfMap)
colnames(Y)[1] <- "Z"

p <- ggplot(Y, aes(x=lon, y=lat)) + theme_bw()
p <- p + theme(plot.title = element_text(size = rel(1.5)))
p <- p + geom_point(aes(colour = Z))
p <- p + scale_colour_gradient(low="pink", 
                               high="black",
                               name = "Z")
p


#### model: RW-BP SOD ####
library(FastKNN)
# matrix of neighbours
nn = matrix(0, N, K) # n x k
for (i in 1:N)
  nn[i,] = k.nearest.neighbors(i, d, k = K)

cluster_nub <- 6:11
cluster_sum <- sum(cluster_nub)
alpha <- 2

C <- list()
center <- c()
for (i in 1:length(cluster_nub)){
  C[[i]] <- kmeans(Y$Z, cluster_nub[i])
  center <- c(center, C[[i]]$centers)
}
Z <- as.matrix(Z)
E <- 1/exp((abs(Z %*% rep(1, cluster_sum) - center))^alpha)
M <- rbind(cbind(t(E), matrix(0, cluster_sum,cluster_sum)),cbind(matrix(0, N,N), E))
W <- matrix(0, N+cluster_sum, N+cluster_sum)
for (i in 1:(N+cluster_sum)){
  for (j in 1:(N+cluster_sum)){
    W[i,j] <- M[i,j]/sum(M[,i])
  }
}

damp_c <- 0.9
S <- matrix(0, N+cluster_sum, N+cluster_sum)
(1-damp_c)*solve(diag(N+cluster_sum)-damp_c*W) %*% 





