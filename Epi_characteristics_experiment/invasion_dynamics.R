library(tidyverse)
library(khroma)
library(ggpubr)

plot_colors <- color("muted")(2)

# contact rate range ------------------------------------------------------

mean.R0s <- data.frame(contact.nums = contact.nums, 
                       H1N1 = H1N1.R0s[1,], 
                       H3N2 = H3N2.R0s[1,])
mean.R0s <- mean.R0s %>%
  pivot_longer(cols=2:3, names_to = "Virus", values_to = "R0")

## calculate 95% confidence intervals
H1N1.CIs <- data.frame(contact.nums = contact.nums,
                  lower = H1N1.R0s[2,], 
                  upper = H1N1.R0s[3,], 
                  Virus = "H1N1")
H3N2.CIs <- data.frame(contact.nums = contact.nums,
                       lower = H3N2.R0s[2,], 
                       upper = H3N2.R0s[3,], 
                       Virus = "H3N2")

R0.CIs <- rbind(H1N1.CIs, H3N2.CIs) %>%
  merge(mean.R0s)

panel_a <- ggplot(R0.CIs, aes(x=contact.nums, y=R0, ymin=lower, ymax=upper, fill=Virus, color=Virus, group=Virus)) +
  geom_point(size=2) +
  geom_errorbar() + 
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  guides(fill="none", colour = guide_legend(override.aes = list(fill=NA,
                                                                linetype=c(0,0)))) +
  scale_y_continuous(limits=c(0, 4), breaks=seq(0, 4, 1)) +
  labs(x="Average number of contacts per day", y=expression(R[0])) +
  theme_light() +
  theme(legend.position = "bottom")

# extinction probability --------------------------------------------------

find.prob.extinction <- function(R0, k){
  ## for a negative binomial offspring distribution
  pgf <- function(s){abs(s-((1+(R0/k)*(1-s))^(-k)))}
  prob <- optimize(pgf, c(0,1))$minimum
  return(prob)
}

## calculate prob of extinction for each subtype at each contact rate
## assuming k=1 and k=large
H1N1.extinction.probs <- matrix(data=NA, nrow=2, ncol=length(contact.nums))
H3N2.extinction.probs <- matrix(data=NA, nrow=2, ncol=length(contact.nums))

for (c in 1:length(contact.nums)){
  H1N1.extinction.probs[1, c] <- find.prob.extinction(H1N1.R0s[1,c], 1)
  H1N1.extinction.probs[2, c] <- find.prob.extinction(H1N1.R0s[1,c], 1e10)
  
  H3N2.extinction.probs[1, c] <- find.prob.extinction(H3N2.R0s[1,c], 1)
  H3N2.extinction.probs[2, c] <- find.prob.extinction(H3N2.R0s[1,c], 1e10)
}

extinction.probs <- data.frame(contact.nums = contact.nums,
                               H1.k1.prob.extinction = H1N1.extinction.probs[1,], 
                               H1.kinf.prob.extinction = H1N1.extinction.probs[2,], 
                               H3.k1.prob.extinction = H3N2.extinction.probs[1,], 
                               H3.kinf.prob.extinction = H3N2.extinction.probs[2,])
extinction.probs <- extinction.probs %>%
  pivot_longer(cols=2:5, names_to="k", values_to="prob.extinction")

panel_b <- ggplot(extinction.probs, aes(x=contact.nums, y=1-prob.extinction, color=k, linetype=k)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[1], plot_colors[1], plot_colors[2], plot_colors[2])) +
  scale_linetype_manual(values=c(2, 1, 2, 1)) +
  guides(color="none", linetype="none") +
  labs(x="Average number of contacts per day", y="Probability of establishment") +
  theme_light()

# stuttering chains -------------------------------------------------------

find.chain.length <- function(R0, k){
  if (R0 < 1){
    mu <- 1 / (1 - R0)
    cov <- sqrt((R0*(1 + (R0/k)))/(1 - R0))
    return(c(mu, cov))
  }
}

## calculate size of stuttering chain at all contact rates where R0 < 1
H1N1.chain.length <- matrix(data=NA, nrow=2, ncol=length(contact.nums), dimnames=list(c("mu", "cov"), c(contact.nums)))
for (c in 1:length(contact.nums)){
  if (H1N1.R0s[1,c] < 1){
    H1N1.chain.length[,c] <- find.chain.length(H1N1.R0s[1,c], H1N1.ks[c])
  }
}

H3N2.chain.length <- matrix(data=NA, nrow=2, ncol=length(contact.nums), dimnames=list(c("mu", "cov"), c(contact.nums)))
for (c in 1:length(contact.nums)){
  if (H3N2.R0s[1,c] < 1){
    H3N2.chain.length[,c] <- find.chain.length(H3N2.R0s[1,c], H3N2.ks[c])
  }
}

chain.lengths <- data.frame(contact.nums = c(c(colnames(H1N1.chain.length[, !colSums(is.na(H1N1.chain.length))])), 
                                             c(colnames(H3N2.chain.length[, !colSums(is.na(H3N2.chain.length))]))),
                            Subtype = c(rep("H1N1", length(H1N1.chain.length[1, !colSums(is.na(H1N1.chain.length))])),
                                        rep("H3N2", length(H3N2.chain.length[1, !colSums(is.na(H3N2.chain.length))]))), 
                            mu = c(H1N1.chain.length[1, !colSums(is.na(H1N1.chain.length))], 
                                   H3N2.chain.length[1, !colSums(is.na(H3N2.chain.length))]), 
                            cov = c(H1N1.chain.length[2, !colSums(is.na(H1N1.chain.length))], 
                                    H3N2.chain.length[2, !colSums(is.na(H3N2.chain.length))]))

panel_c <- ggplot(chain.lengths, aes(x=as.numeric(contact.nums), y=mu, color=Subtype, group=Subtype)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), name="Virus") + 
  scale_x_continuous(breaks=seq(5, 25, 5), limits=c(5, 25)) +
  scale_y_continuous(breaks=seq(0, 40, 5), limits=c(0, 37)) +
  guides(color="none") +
  labs(x="Average number of contacts per day", y="Average length of stuttering chain") +
  theme_light() 

# intrinsic growth rate ---------------------------------------------------

find.growth.rate.exp <- function(R, Tc){
  exponential.r <- (R-1) / Tc
  return(exponential.r)
}

find.growth.rate.delta <- function(R, Tc){
  delta.r <- log(R) / Tc
  return(delta.r)
}

## using n=5
find.growth.rate.gamma <- function(R, Tc){
  gamma.r <- (5/Tc)*R^(1/5) - (5/Tc)
  return(gamma.r)
}

## assume intermediate gen interval 
H1N1.Tc <- 3.6
H3N2.Tc <- 3.6

H1N1.growth.rates <- data.frame(Virus = rep("H1N1", length(H1N1.R0s[1,])), 
                                exponential = find.growth.rate.exp(H1N1.R0s[1,], H1N1.Tc), 
                                delta = find.growth.rate.delta(H1N1.R0s[1,], H1N1.Tc), 
                                gamma = find.growth.rate.gamma(H1N1.R0s[1,], H1N1.Tc),
                                contact.rate = contact.nums)
H3N2.growth.rates <- data.frame(Virus = rep("H3N2", length(H3N2.R0s[1,])), 
                                exponential = find.growth.rate.exp(H3N2.R0s[1,], H3N2.Tc), 
                                delta = find.growth.rate.delta(H3N2.R0s[1,], H3N2.Tc), 
                                gamma = find.growth.rate.gamma(H3N2.R0s[1,], H3N2.Tc),
                                contact.rate = contact.nums)

combined.growth.rates <- rbind(H1N1.growth.rates, H3N2.growth.rates)

combined.growth.rates <- combined.growth.rates %>%
  pivot_longer(cols=2:4, values_to="growth.rate", names_to="dist")

panel_d <- ggplot(combined.growth.rates, aes(x=contact.rate, y=growth.rate, group=interaction(Virus, dist), color=Virus, linetype=dist)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_linetype_manual(values = c(2, 3, 1), name="Distribution") +
  guides(color="none") +
  labs(x="Average number of contacts per day", y="Intrinsic growth rate") +
  theme_light() +
  theme(legend.position= "inside", legend.position.inside = c(0.25, 0.8),
        legend.key.width = unit(4, "line"), legend.background = element_rect(fill = "white", color = "black"))

# plot -------------------------------------------------------------------

top <- ggarrange(panel_a, panel_b, ncol=2, align="h", common.legend = T, legend = "top", labels=c("A", "B"))
bottom <- ggarrange(panel_c, panel_d, ncol=2, align="h", labels=c("C", "D"))
ggarrange(top, bottom, nrow=2)
