library(tidyverse)
library(khroma)
library(ggpubr)

# extinction probability --------------------------------------------------

H1N1.R0 <- 2.162737
H1N1.k <- 174.6621
H3N2.R0 <- 0.4631673
H3N2.k <- 4.721216

prob.extinction <- function(R0, k){
  ## for a negative binomial offspring distribution
  pgf <- function(s){abs(s-((1+(R0/k)*(1-s))^(-k)))}
  prob <- optimize(pgf, c(0,1))$minimum
  return(prob)
}

R0 <- 1.63
k <- 1 #2359926

x <- seq(0, 1, 0.001)
y <- (1 + (R0/k)*(1-x))^-k
plot(x, y, ylim=c(0, 1)) 
abline(a=0, b=1)

prob.extinction(R0, k)

extinction.probs <- data.frame(Subtype = c("H1N1", "H3N2"), 
                               prob = c(prob.extinction(H1N1.R0, H1N1.k), prob.extinction(H3N2.R0, H3N2.k)))

#prob.extinction(H1N1.R0, H1N1.k)
#prob.extinction(H3N2.R0, H3N2.k)

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

length.chains <- data.frame(Subtype = "H3N2", 
                            length = chain.length(H3N2.R0, H3N2.k)[[1]], 
                            cov = chain.length(H3N2.R0, H3N2.k)[[2]])

# intrinsic growth rate ---------------------------------------------------

H1N1.shape <- 2.9592369
H1N1.rate <- 0.7198616
#H1N1.G <- 4.11
#H1N1.sigma <- sqrt(2.9592369/(0.7198616^2))
H3N2.shape <- 2.6401310
H3N2.rate <- 0.6038136
#H3N2.G <- 4.37
#H3N2.sigma <- sqrt(2.6401310/(0.6038136^2))

growth.rate <- function(R, G, sigma){
  r <- - (-G + sqrt(G^2 - 2*sigma^2*log(R))) / sigma^2
  return(r)
}

rates <- data.frame(Subtype = c("H1N1", "H3N2"), 
                    r = c(growth.rate(H1N1.R0, H1N1.G, H1N1.sigma), growth.rate(H3N2.R0, H3N2.G, H3N2.sigma)))

#growth.rate(H1N1.R0, H1N1.G, H1N1.sigma)
#growth.rate(H3N2.R0, H3N2.G, H3N2.sigma)

rates <- data.frame(Subtype = c("H1N1", "H3N2"), 
                    r = c(H1N1.rate*(H1N1.R0^(1/H1N1.shape) - 1), H3N2.rate*(H3N2.R0^(1/H3N2.shape) - 1)))

H1N1.rate*(H1N1.R0^(1/H1N1.shape) - 1)
H3N2.rate*(H3N2.R0^(1/H3N2.shape) - 1)

# plots -------------------------------------------------------------------

plot_colors <- color("muted")(2)

panel_a <- ggplot(extinction.probs, aes(x=Subtype, y=1-prob, color=Subtype, fill=Subtype)) +
  geom_col(width=0.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  guides(color="none", fill="none") +
  labs(x=NULL, y="Probability of establishment") +
  theme_light()

ggplot(length.chains, aes(x=Subtype, y=length)) +
  geom_point(size=4, color=plot_colors[[2]]) +
  geom_segment(aes(x=Subtype, xend=Subtype, y=0, yend=length), linewidth = 1.5, color=plot_colors[[2]]) +
  labs(x="Subtype", y="Average length of stuttering chain") +
  theme_light()

panel_b <- ggplot(rates, aes(x=Subtype, y=r, color=Subtype, fill=Subtype)) +
  geom_col(width=0.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  guides(color="none", fill="none") +
  ylim(-0.25, 0.25) +
  labs(x=NULL, y="Intrinsic growth rate") +
  theme_light()

ggarrange(panel_a, panel_b, ncol=2, align="h", labels=c("A", "B"))

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
