library(khroma)
library(tidyverse)
library(ggpubr)
library(DescTools)

# panel A -----------------------------------------------------------------

## choose representative example donor/contact pair; here from H1N1 experiment
ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "dpi", "nw_titer", "donor_dose")

pair_5444_5445 <- ferrets %>%
  filter(Ferret_ID == "5444" | Ferret_ID == "5445") %>%
  mutate(dpe = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

panel_a <- ggplot(pair_5444_5445, aes(x=dpe, y=nw_titer, color=DI_RC, shape=LOD_shape)) +
  geom_point(size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1.5) +
  scale_color_manual(labels = c("Donor", "Recipient"), values = c("black", "red")) +
  scale_shape_manual(values=c(16, 21)) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme_light() +
  theme(legend.position = "top") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  labs(title = NULL, x = NULL, y = expression(paste("Viral titer (", log[10], TCID[50], ")")), color=NULL)

# panel B -----------------------------------------------------------------

## transmission probabilities at varying s values

load("panel_3b_data.RData")

s_0.05 <- panel_3b[[1]][["5444"]]
s_0.05 <- s_0.05 %>%
  select(dpe, pr_contact_pos_log) %>%
  mutate(s_val = "0.05")
s_0.1 <- panel_3b[[2]][["5444"]]
s_0.1 <- s_0.1 %>%
  select(dpe, pr_contact_pos_log) %>%
  mutate(s_val = "0.1")
s_0.2 <- panel_3b[[3]][["5444"]]
s_0.2 <- s_0.2 %>%
  select(dpe, pr_contact_pos_log) %>%
  mutate(s_val = "0.2")
combined_df <- rbind(s_0.05, s_0.1, s_0.2)
## set initial probability of transmission to 0
combined_df$pr_contact_pos_log <- replace_na(combined_df$pr_contact_pos_log, 0)

## binary contact data
partner_df <- pair_5444_5445 %>%
  filter(Ferret_ID == "5444") %>%
  ## the contact animal is or has previously been positive every time point after 0 dpe
  mutate(pos_neg = c(0, 1, 1, 1, 1, 1)) %>%
  mutate(point = if_else(pos_neg < 0.5, "below", "above"))

plot_colors <- color("muted", reverse=T)(3)

## need to make multiple line plot (split lines)

panel_b <- ggplot(combined_df, aes(x=dpe, y=pr_contact_pos_log, color=s_val)) +
  geom_line(linewidth=1.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]], plot_colors[[3]])) + 
  geom_point(data=partner_df, aes(x=dpe, y=pos_neg, shape=point), color="red", size=2) +
  scale_shape_manual(values = c(16, 1)) +
  guides(shape="none") +
  labs(x="Days post exposure", y="Cumulative transmission probability", shape=NULL, color="s value") +
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 10, 2))

combo_test_df$prob <- replace_na(combo_test_df$prob, 0)

multi_line_s <- ggplot(combo_test_df, aes(x=dpe, y=prob, group=interaction(before_cut, s_val), color=s_val, linetype=before_cut)) +
  geom_line(linewidth=1.5) +
  scale_linetype_manual(values=c("solid", "twodash")) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]], plot_colors[[3]])) +
  labs(x="Days post exposure", y="Cumulative transmission probability", shape=NULL, color="s value") +
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 10, 2)) +
  geom_point(data=partnerdf, aes(x=dpch, y=pos_neg, shape=Point), inherit.aes = F, color="red", size=2) +
  scale_shape_manual(values = c(16, 1)) +
  guides(shape="none", linetype="none")


# panel C -----------------------------------------------------------------

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

panel_c <- ggplot(combo_MLE_trace, aes(x=s, y=likelihood, color=virus)) +
  geom_line(linewidth=1.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  labs(x="s", y="Log likelihood", color="Virus") +
  geom_vline(xintercept=MLE_H1N1, color=plot_colors[1], linewidth=1, linetype = 2) +
  geom_vline(xintercept=MLE_H3N2, color=plot_colors[2], linewidth=1, linetype = 2) +
  annotate("rect", xmin=CIs_H1N1[1], xmax=CIs_H1N1[2], ymin=min(combo_MLE_trace$likelihood), ymax=max(combo_MLE_trace$likelihood), alpha=0.2, color=plot_colors[1], fill=plot_colors[1]) +
  annotate("rect", xmin=CIs_H3N2[1], xmax=CIs_H3N2[2], ymin=min(combo_MLE_trace$likelihood), ymax=max(combo_MLE_trace$likelihood), alpha=0.2, color=plot_colors[2], fill=plot_colors[2]) +
  theme_light() +
  theme(legend.position = "top")


# panel D -----------------------------------------------------------------

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
    prob <- calculate_pr_constant(MLE_H1N1, slice$log_VL)
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
    prob <- calculate_pr_constant(MLE_H3N2, slice$log_VL)
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
    prob_upper <- calculate_pr_constant(CIs_H1N1[2], slice$log_VL)
    prob_lower <- calculate_pr_constant(CIs_H1N1[1], slice$log_VL)
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
    prob_upper <- calculate_pr_constant(CIs_H3N2[2], slice$log_VL)
    prob_lower <- calculate_pr_constant(CIs_H3N2[1], slice$log_VL)
  } else {
    prob_upper <- 0
    prob_lower <- 0
  }
  H3N2_ribbon$prob_upper[row] <- prob_upper
  H3N2_ribbon$prob_lower[row] <- prob_lower
}

plot_colors <- color("muted")(2)

panel_d <- ggplot(combo_VL_pr, aes(x=log_VL, y=pr, color=virus)) +
  geom_line(linewidth=1.5) +
  geom_ribbon(data=H1N1_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[1]], alpha = 0.2, inherit.aes = F) +
  geom_ribbon(data=H3N2_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[2]], alpha = 0.2, inherit.aes = F) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  labs(x=expression(paste("Viral titer (", log[10], TCID[50], ")")), y="Probability of transmission", color="Virus") +
  theme_light() +
  guides(color="none")


# full figure -------------------------------------------------------------

## combine

ggarrange(panel_a, panel_c, panel_b, panel_d, nrow=2, ncol=2, labels = c("A", "C", "B", "D"), widths = c(1.4, 1), align = "v")
