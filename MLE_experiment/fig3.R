library(khroma)
library(tidyverse)
library(ggpubr)
library(DescTools)

# panel A -----------------------------------------------------------------

## results from MLE for H1N1
load("MLE_experiment/H1N1_MLE_trace.RData")
MLE_H1N1 <- 0.111
CIs_H1N1 <- c(0.068, 0.172)

## results from MLE for H3N2
load("MLE_experiment/H3N2_MLE_trace.RData")
MLE_H3N2 <- 0.044
CIs_H3N2 <- c(0.021, 0.079)

combo_MLE_trace <- merge(H1N1_MLE_trace, H3N2_MLE_trace, by="s", suffixes = c(".H1N1", ".H3N2"))
combo_MLE_trace <- combo_MLE_trace %>%
  pivot_longer(-s, names_to="virus", names_prefix = "prob.", values_to = "likelihood")
## remove probabilities at s=0
combo_MLE_trace <- combo_MLE_trace[3:length(combo_MLE_trace$s),]

plot_colors <- color("muted")(2)

panel_a <- ggplot(combo_MLE_trace, aes(x=s, y=likelihood, color=virus)) +
  geom_line(linewidth=1.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  labs(x="s", y="Log likelihood", color="Subtype") +
  geom_vline(xintercept=MLE_H1N1, color=plot_colors[1], linewidth=1, linetype = 2) +
  geom_vline(xintercept=MLE_H3N2, color=plot_colors[2], linewidth=1, linetype = 2) +
  annotate("rect", xmin=CIs_H1N1[1], xmax=CIs_H1N1[2], ymin=min(combo_MLE_trace$likelihood), ymax=max(combo_MLE_trace$likelihood), alpha=0.2, color=plot_colors[1], fill=plot_colors[1]) +
  annotate("rect", xmin=CIs_H3N2[1], xmax=CIs_H3N2[2], ymin=min(combo_MLE_trace$likelihood), ymax=max(combo_MLE_trace$likelihood), alpha=0.2, color=plot_colors[2], fill=plot_colors[2]) +
  theme_light() +
  theme(legend.position = "top")

# panel B -----------------------------------------------------------------

## for Pr+ive

calculate_pr_constant <- function(s, VL){
  lambda <- s * VL
  ## integrate from 0 to 1
  integral <- AUC(x=c(0,1), y=c(lambda, lambda), method = "trapezoid")
  prob <- 1 - exp(-integral)
  return (prob)
}

VL_probs_H1N1 <- data.frame(log_VL = seq(0, 10, 0.1), 
                            prob = rep(0, 101))

for (row in 1:nrow(VL_probs_H1N1)){
  slice <- VL_probs_H1N1[row,]
  if (slice$log_VL > 0.5){
    prob <- 1 - exp(-slice$log_VL * MLE_H1N1)
  } else {
    prob <- 0
  }
  VL_probs_H1N1$prob[row] <- prob
}

VL_probs_H3N2 <- data.frame(log_VL = seq(0, 10, 0.1), 
                            prob = rep(0, 101))
for (row in 1:nrow(VL_probs_H3N2)){
  slice <- VL_probs_H3N2[row,]
  if (slice$log_VL > 0.5){
    prob <- 1 - exp(-slice$log_VL * MLE_H3N2)
  } else {
    prob <- 0
  }
  VL_probs_H3N2$prob[row] <- prob
}

combo_VL_pr <- merge(VL_probs_H1N1, VL_probs_H3N2, by="log_VL", suffixes = c(".H1N1", ".H3N2"))
combo_VL_pr <- combo_VL_pr %>%
  pivot_longer(-log_VL, names_to="virus", names_prefix = "prob.", values_to = "pr")

## now, calculate CIs for ribbon
H1N1_ribbon <- data.frame(log_VL = seq(0, 10, 0.1), 
                          prob_upper = rep(0, 101), 
                          prob_lower = rep(0, 101))
for (row in 1:nrow(H1N1_ribbon)){
  slice <- H1N1_ribbon[row,]
  if (slice$log_VL > 0.5){
    prob_upper <- 1 - exp(-slice$log_VL * CIs_H1N1[2])
    prob_lower <- 1 - exp(-slice$log_VL * CIs_H1N1[1])
  } else {
    prob_upper <- 0
    prob_lower <- 0
  }
  H1N1_ribbon$prob_upper[row] <- prob_upper
  H1N1_ribbon$prob_lower[row] <- prob_lower
}

H3N2_ribbon <- data.frame(log_VL = seq(0, 10, 0.1), 
                          prob_upper = rep(0, 101), 
                          prob_lower = rep(0, 101))
for (row in 1:nrow(H3N2_ribbon)){
  slice <- H3N2_ribbon[row,]
  if (slice$log_VL > 0.5){
    prob_upper <- 1 - exp(-slice$log_VL * CIs_H3N2[2])
    prob_lower <- 1 - exp(-slice$log_VL * CIs_H3N2[1])
  } else {
    prob_upper <- 0
    prob_lower <- 0
  }
  H3N2_ribbon$prob_upper[row] <- prob_upper
  H3N2_ribbon$prob_lower[row] <- prob_lower
}

plot_colors <- color("muted")(2)

panel_b <- ggplot(combo_VL_pr, aes(x=log_VL, y=pr, color=virus)) +
  geom_line(linewidth=1.5) +
  geom_ribbon(data=H1N1_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[1]], alpha = 0.2, inherit.aes = F) +
  geom_ribbon(data=H3N2_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[2]], alpha = 0.2, inherit.aes = F) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  labs(x=expression(paste("Viral titer (", log[10], TCID[50], ")")), y="Probability of transmission", color="Virus") +
  theme_light() +
  guides(color="none")

# full figure -------------------------------------------------------------

## combine

ggarrange(panel_a, panel_b, ncol=2, labels = c("A", "B"), common.legend = T, legend = "bottom")
