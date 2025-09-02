# em algorithm: https://teng-gao.github.io/blog/2022/ems/
library(dplyr)

get.phi.k = function(k, lambda.t, pi.t, x.i) { 
  l = sapply(1:K, function(k){
    lambda.t[k] * exp(-lambda.t[k] * x.i) * pi.t[k]
  })
  
  phi.k = l[k]/sum(l) 
  return(phi.k)
}

update.lambda.k = function(k, phi.k, x){
  lamda.k = sum(phi.k)/ sum(phi.k * x)
  return(lamda.k)
}

update.pi.k = function(k, phi.k){
  pi.k = sum(phi.k) / N
  return(pi.k)
}


run.em <- function(x, max.iter = 10){
  
  # initialize
  pi <- list(); lambda <- list()
  pi[[1]] <- rep(1/K, K)
  lambda[[1]] <- 1:K
  
  for (t in 1:(max.iter-1)){
    
    # E-step
    phi = lapply(1:K, function(k){
      sapply(1:N, function(i){
          get.phi.k(k=k, lambda.t = lambda[[t]], pi.t = pi[[t]], x.i=x)
      })
    })
    
    
    # M-step
    lambda[[t+1]] <- sapply(1:K, function(k){
      update.lambda.k(k=k, phi.k = phi[[k]], x = x)
    })
    
    pi[[t+1]] <- sapply(1:K, function(k){
      update.pi.k(k=k, phi.k=phi[[k]])
    })
  }
  
  return(list(lambda = lambda, pi = pi))
  
}


# test
K <- 3
N <- 1000
set.seed(2)

x <- c(
    rgamma(n = 500, shape = 1, rate = 1),
    rgamma(n = 300, shape = 1, rate = 10),
    rgamma(n = 200, shape = 1, rate = 100)
) %>% sort

n.iter <- 100
res <- run.em(x = x, max.iter = n.iter)