library(tidyverse)
library(khroma)
library(ggpubr)

plot_colors <- color("muted")(2)


# LINEAR ------------------------------------------------------------------

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
  labs(x="Number of contacts per day", y=expression(R[0])) +
  theme_light() +
  theme(legend.position = "bottom")


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
  labs(x="Number of contacts per day", y="Probability of establishment") +
  theme_light()


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
  scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0, 25)) +
  guides(color="none") +
  labs(x="Number of contacts per day", y="Average length of stuttering chain") +
  theme_light() 


## params chosen from mean gen time/cov scatter
find.growth.rate.gamma <- function(R, alpha=10.822, lambda=3.7911){
  gamma.r <- lambda * (R^(1/alpha) - 1)
  return(gamma.r)
}

H1N1.growth.rate <- find.growth.rate.gamma(H1N1.R0s[1,])
H3N2.growth.rate <- find.growth.rate.gamma(H3N2.R0s[1,])
combined.growth.rates <- data.frame(Virus = c(rep("H1N1", length(contact.nums)), rep("H3N2", length(contact.nums))),
                                    contact.rate = c(contact.nums, contact.nums),
                                    rate = c(H1N1.growth.rate, H3N2.growth.rate))

panel_d <- ggplot(combined.growth.rates, aes(x=contact.rate, y=rate, color=Virus)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  guides(color="none") +
  labs(x="Number of contacts per day", y=expression(paste("Intrinsic growth rate, r (", days^-1, ")"))) +
  theme_light()

top <- ggarrange(panel_a, panel_b, ncol=2, align="h", common.legend = T, legend = "top", labels=c("A", "B"))
bottom <- ggarrange(panel_c, panel_d, ncol=2, align="h", labels=c("C", "D"))
ggarrange(top, bottom, nrow=2)


# THRESHOLD ---------------------------------------------------------------

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
  labs(x="Number of contacts per day", y=expression(R[0])) +
  theme_light() +
  theme(legend.position = "bottom")


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
  labs(x="Number of contacts per day", y="Probability of establishment") +
  theme_light()


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
  scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0, 25)) +
  guides(color="none") +
  labs(x="Number of contacts per day", y="Average length of stuttering chain") +
  theme_light() 


## params chosen from mean gen time/cov scatter
find.growth.rate.gamma <- function(R, alpha=4.53, lambda=1.27){
  gamma.r <- lambda * (R^(1/alpha) - 1)
  return(gamma.r)
}

H1N1.growth.rate <- find.growth.rate.gamma(H1N1.R0s[1,])
H3N2.growth.rate <- find.growth.rate.gamma(H3N2.R0s[1,])
combined.growth.rates <- data.frame(Virus = c(rep("H1N1", length(contact.nums)), rep("H3N2", length(contact.nums))),
                                    contact.rate = c(contact.nums, contact.nums),
                                    rate = c(H1N1.growth.rate, H3N2.growth.rate))

panel_d <- ggplot(combined.growth.rates, aes(x=contact.rate, y=rate, color=Virus)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  guides(color="none") +
  labs(x="Number of contacts per day", y=expression(paste("Intrinsic growth rate, r (", days^-1, ")"))) +
  theme_light()


top <- ggarrange(panel_a, panel_b, ncol=2, align="h", common.legend = T, legend = "top", labels=c("A", "B"))
bottom <- ggarrange(panel_c, panel_d, ncol=2, align="h", labels=c("C", "D"))
ggarrange(top, bottom, nrow=2)

# HILL --------------------------------------------------------------------

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
  labs(x="Number of contacts per day", y=expression(R[0])) +
  theme_light() +
  theme(legend.position = "bottom")


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
  labs(x="Number of contacts per day", y="Probability of establishment") +
  theme_light()


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
  scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0, 25)) +
  guides(color="none") +
  labs(x="Number of contacts per day", y="Average length of stuttering chain") +
  theme_light() 


## params chosen from mean gen time/cov scatter
find.growth.rate.gamma <- function(R, alpha=2.579, lambda=0.738){
  gamma.r <- lambda * (R^(1/alpha) - 1)
  return(gamma.r)
}

H1N1.growth.rate <- find.growth.rate.gamma(H1N1.R0s[1,])
H3N2.growth.rate <- find.growth.rate.gamma(H3N2.R0s[1,])
combined.growth.rates <- data.frame(Virus = c(rep("H1N1", length(contact.nums)), rep("H3N2", length(contact.nums))),
                                    contact.rate = c(contact.nums, contact.nums),
                                    rate = c(H1N1.growth.rate, H3N2.growth.rate))

panel_d <- ggplot(combined.growth.rates, aes(x=contact.rate, y=rate, color=Virus)) +
  geom_point(size=3) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  guides(color="none") +
  labs(x="Number of contacts per day", y=expression(paste("Intrinsic growth rate, r (", days^-1, ")"))) +
  theme_light()


top <- ggarrange(panel_a, panel_b, ncol=2, align="h", common.legend = T, legend = "top", labels=c("A", "B"))
bottom <- ggarrange(panel_c, panel_d, ncol=2, align="h", labels=c("C", "D"))
ggarrange(top, bottom, nrow=2)
