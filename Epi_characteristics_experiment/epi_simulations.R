library(tidyverse)
library(khroma)
library(ggpubr)

plot_colors <- color("muted")(2)

# contact rate range ------------------------------------------------------

mean.R0s <- data.frame(contact.nums = contact.nums, 
                       H1N1 = H1N1.R0s, 
                       H3N2 = H3N2.R0s)
mean.R0s <- mean.R0s %>%
  pivot_longer(cols=2:3, names_to = "Subtype", values_to = "R0")

panel_a <- ggplot(mean.R0s, aes(x=contact.nums, y=R0, color=Subtype, fill=Subtype)) +
  geom_col(position = "dodge") +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  guides(color="none", fill="none") +
  labs(x="Number of contacts per day", y=expression(R[0])) +
  theme_light()

# extinction probability --------------------------------------------------

find.prob.extinction <- function(R0, k){
  ## for a negative binomial offspring distribution
  pgf <- function(s){abs(s-((1+(R0/k)*(1-s))^(-k)))}
  prob <- optimize(pgf, c(0,1))$minimum
  return(prob)
}

## assume very large k for now??
k <- 1000

## calculate prob of extinction for each subtype at each contact rate
H1N1.extinction.probs <- c()
H3N2.extinction.probs <- c()

for (c in 1:length(contact.nums)){
  H1.R <- H1N1.R0s[c]
  H1.k <- H1N1.ks[c]
  H1N1.extinction.probs <- append(H1N1.extinction.probs, find.prob.extinction(H1.R, H1.k))
  H3.R <- H3N2.R0s[c]
  H3.k <- H3N2.ks[c]
  H3N2.extinction.probs <- append(H3N2.extinction.probs, find.prob.extinction(H3.R, H3.k))
}

extinction.probs <- data.frame(contact.nums = contact.nums,
                               H1N1 = H1N1.extinction.probs, 
                               H3N2 = H3N2.extinction.probs)
extinction.probs <- extinction.probs %>%
  pivot_longer(cols=2:3, names_to="Subtype", values_to="prob.extinction")

panel_b <- ggplot(extinction.probs, aes(x=contact.nums, y=1-prob.extinction, color=Subtype, fill=Subtype)) +
  geom_col(position="dodge") +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  guides(color="none", fill="none") +
  labs(x="Number of contacts", y="Probability of etablishment") +
  theme_light()

# intrinsic growth rate ---------------------------------------------------

find.growth.rate.exp <- function(R, G){
  exponential.r <- (R-1) / G
  return(exponential.r)
}

find.growth.rate.delta <- function(R, G){
  delta.r <- log(R) / G
  return(delta.r)
}

H1N1.growth.rates <- data.frame(Subtype = rep("H1N1", length(H1N1.R0s)), 
                                exponential.r = find.growth.rate.exp(H1N1.R0s, H1N1.Tcs), 
                                delta.r = find.growth.rate.delta(H1N1.R0s, H1N1.Tcs), 
                                contact.rate = contact.nums)
H3N2.growth.rates <- data.frame(Subtype = rep("H3N2", length(H3N2.R0s)), 
                                exponential.r = find.growth.rate.exp(H3N2.R0s, H3N2.Tcs), 
                                delta.r = find.growth.rate.delta(H3N2.R0s, H3N2.Tcs), 
                                contact.rate = contact.nums)

combined.growth.rates <- rbind(H1N1.growth.rates, H3N2.growth.rates)

panel_c <- ggplot(combined.growth.rates, aes(x=contact.rate, y=exponential.r, group=Subtype, color=Subtype)) +
  geom_linerange(aes(ymin=delta.r, ymax=exponential.r), linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  labs(x="Number of contacts", y="Intrinsic growth rate") +
  theme_light()

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

find.chain.length <- function(R0, k, j){
  prob.chain.size
}

length.chains <- data.frame(Subtype = "H3N2", 
                            length = chain.length(H3N2.R0, H3N2.k)[[1]], 
                            cov = chain.length(H3N2.R0, H3N2.k)[[2]])

# plots -------------------------------------------------------------------

ggplot(length.chains, aes(x=Subtype, y=length)) +
  geom_point(size=4, color=plot_colors[[2]]) +
  geom_segment(aes(x=Subtype, xend=Subtype, y=0, yend=length), linewidth = 1.5, color=plot_colors[[2]]) +
  labs(x="Subtype", y="Average length of stuttering chain") +
  theme_light()


ggarrange(panel_a, panel_b, panel_c, nrow=2, ncol=2)

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
