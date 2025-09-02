# em algorithm: https://teng-gao.github.io/blog/2022/ems/
# install.packages("rstudioapi")
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)

setwd("/Users/adamkurth/Documents/vscode/code/understanding-series/understanding-mcmc")

# E-step: calculate posterior probability for single data point
get.phi.k = function(k, lambda.t, pi.t, x.i) { 
  
  # calculate likelihood contribution from each component
  l = sapply(1:K, function(k){
    lambda.t[k] * exp(-lambda.t[k] * x.i) * pi.t[k]
  })
  
  # avoid NaN
  if (sum(l) == 0) return(1/K)
  
  phi.k = l[k]/sum(l) # return for component k
  return(phi.k)
}


# M-Step: Update the rate parameter lambda
update.lambda.k = function(k, phi.k, x){
  
  # add small epsilon prevent division by 0
  lamda.k = sum(phi.k)/ (sum(phi.k * x) + 1e-10)
  return(lamda.k)
}

# M-step: update mixing proportion pi
update.pi.k = function(k, phi.k){
  pi.k = sum(phi.k) / N
  return(pi.k)
}

# main
run.em <- function(x, max.iter = 10){
  
  # initialize
  pi <- list(); lambda <- list()
  pi[[1]] <- rep(1/K, K)
  lambda[[1]] <- 1:K
  
  # data frame tidy format
  results.df <- data.frame(
    iter = 1, 
    k = 1:K, 
    lambda = lambda[[1]], 
    pi = pi[[1]]
  )
  
  for (t in 1:(max.iter-1)){
    
    # --- E-step ---
    # calculate posterior posterior for all data points for each component
    
    phi <- lapply(1:K, function(k){
      sapply(1:N, function(i){
          get.phi.k(k=k, lambda.t = lambda[[t]], pi.t = pi[[t]], x.i=x[i])
      })
    })
    
    
    # --- M-step ---
    # Update lambda and pi for each component
    
    lambda.new <- sapply(1:K, function(k){
      update.lambda.k(k = k, phi.k = phi[[k]], x = x)
    })
    
    pi.new <- sapply(1:K, function(k){
      update.pi.k(k = k, phi.k = phi[[k]])
    })
    
    # store
    lambda[[t+1]] <- lambda.new
    pi[[t+1]] <- pi.new 
    
    # append
    results.df <- rbind(results.df, data.frame(
      iter = t+1, 
      k = 1:K, 
      lambda = lambda.new,
      pi = pi.new
    ))
  }
  
  return(results.df)
  
}


# data generaiton / execution
K <- 3
N <- 1000
set.seed(2)

# mixture of 3 exp distributions
x <- c(
    rgamma(n = 500, shape = 1, rate = 1),
    rgamma(n = 300, shape = 1, rate = 10),
    rgamma(n = 200, shape = 1, rate = 100)
) %>% sort

n.iter <- 100
em.res <- run.em(x = x, max.iter = n.iter)







# --- Visualization with ggplot2 and gganimate ---

# 1. Prepare data for the density curves
density_data <- em.res %>%
  group_by(iter) %>%
  summarise(
    # Create a grid of x values for the plot
    x_vals = list(seq(0, max(x), length.out = 500)),
    # Calculate the mixture density for each x value
    y_vals = list(
      rowSums(
        sapply(1:K, function(k_idx) {
          pi[k_idx] * dexp(x_vals[[1]], rate = lambda[k_idx])
        })
      )
    )
  ) %>%
  tidyr::unnest(c(x_vals, y_vals))

# 2. Prepare data for the text labels
text_data <- em.res %>%
  group_by(iter) %>%
  summarise(
    label = paste(
      # Use plotmath syntax: '==' for equals, '[]' for subscripts, '*' for spacing
      sprintf("lambda[%d] == %.2f*','~pi[%d] == %.2f", k, lambda, k, pi), 
      collapse = "\n"
    )
  )

# Create the animation
p <- ggplot() +
  geom_histogram(data = data.frame(x = x), aes(x = x, y = ..density..),
                 bins = 50, fill = "lightblue", color = "grey80") +
  geom_line(data = density_data, aes(x = x_vals, y = y_vals), color = "red", size = 1.2) +
  
  # 2. Add `parse = TRUE` to geom_text to render the expressions
  geom_text(data = text_data, aes(x = Inf, y = Inf, label = label),
            parse = TRUE, # This is the crucial change!
            hjust = 1.05, vjust = 1.1, lineheight = 1.1, size = 4) +
  
  labs(
    title = 'EM Algorithm for Mixture of Exponentials',
    subtitle = 'Iteration: {closest_state}',
    x = 'x',
    y = 'Density'
  ) +
  theme_minimal(base_size = 14) +
  ylim(0, 10) +
  transition_states(iter, transition_length = 2, state_length = 1)

anim <- animate(p, nframes = n.iter, fps = 2, width = 800, height = 500)
anim_save("em.gif", animation = anim)
anim

