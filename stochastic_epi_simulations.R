library(tidyverse)
library(khroma)

load("H1N1_distr.RData")
load("H3N2_distr.RData")

# extinction probability --------------------------------------------------

prob.extinction <- function(R0, k){
  ## for a negative binomial offspring distribution
  pgf <- function(s){abs(s-((1+(R0/k)*(1-s))^(-k)))}
  prob <- optimize(pgf, c(0,1))$minimum
  return(prob)
}

R0 <- 1.63
k <- 2359926

prob.extinction(R0, k)

## code from Madison for calculating the probability of extinction given an R0 and k

extinction_prob <- function(R0, k){ 
  if (is.infinite(k)){
    #pgf = exp(-R0*(1-s)), constant R0
    getProb <- function(x){abs(x-exp(-R0*(1-x)))} 
  }
  else {
    #pgf = (1+(R0/k)*(1-s))^(-k)
    getProb <- function(x){abs(x-(1+(R0/k)*(1-x))^(-k))}
    #(1/R0)*integral(fun,0,x)-0.2)); 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2
  }
  extinction.prob <- optimize(getProb, c(0, 1))$minimum
  return(extinction.prob)
}

extinction_prob(R0, k)

# stuttering chains -------------------------------------------------------

chain.length <- function(R0, k){
  if (R0 > 1){
    print("R0 is greater than 1")
  } else {
    mu <- 1 / (1 - R0)
    cov <- sqrt((R0*(1 + (R0/k)))/(1 - R0))
    return(c("mu" = mu, "cov" = cov))
  }
}

H3N2_chains <- chain.length(H3N2_mean_R0, H3N2_mean_k)

# branching process simulations ---------------------------------------------------------

n.simulations <- 100
generations <- 25

H1N1_simulations <- matrix(0, nrow=generations, ncol=n.simulations)
H1N1.extinct <- 0

## using neg binomial distribution
for (n in 1:n.simulations){
  pop_size <- c(1) # start with one person
  for (gen in 2:generations){
    Znplus1 <- rnbinom(tail(pop_size, n=1), size=H1N1_k, mu=H1N1_mean_R0)
    pop_size <- append(pop_size, sum(Znplus1)) ## add the sum of the draws for total size
  }
  if (tail(pop_size, n=1) == 0){
    H1N1.extinct <- H1N1.extinct + 1
  }
  H1N1_simulations[1:length(pop_size), n] <- pop_size
}

H1N1_plot_sims <- as.data.frame(H1N1_simulations) %>%
  mutate(time = 1:generations)
H1N1_plot_sims <- H1N1_plot_sims %>%
  pivot_longer(cols=1:n.simulations, names_to="simulation", values_to="pop_size")

ggplot(H1N1_plot_sims, aes(x=time, y=pop_size, color=simulation)) +
  geom_line() +
  guides(colour = "none") + 
  theme_light() +
  ggtitle(paste("simulations for nbinom, mean=", H1N1_mean_R0, "shape=", H1N1_k))


H3N2_simulations <- matrix(0, nrow=generations, ncol=n.simulations)
H3N2.extinct <- 0

## using neg binomial distribution
for (n in 1:n.simulations){
  pop_size <- c(1) # start with one person
  for (gen in 2:generations){
    Znplus1 <- rnbinom(tail(pop_size, n=1), size=H3N2_k, mu=H3N2_mean_R0)
    pop_size <- append(pop_size, sum(Znplus1)) ## add the sum of the draws for total size
  }
  if (tail(pop_size, n=1) == 0){
    H3N2.extinct <- H3N2.extinct + 1
  }
  H3N2_simulations[1:length(pop_size), n] <- pop_size
}

H3N2_plot_sims <- as.data.frame(H3N2_simulations) %>%
  mutate(time = 1:generations)
H3N2_plot_sims <- H3N2_plot_sims %>%
  pivot_longer(cols=1:n.simulations, names_to="simulation", values_to="pop_size")

ggplot(H3N2_plot_sims, aes(x=time, y=pop_size, color=simulation)) +
  geom_line() +
  guides(colour = "none") + 
  theme_light() +
  ggtitle(paste("simulations for nbinom, mean=", H3N2_mean_R0, "shape=", H3N2_k))
