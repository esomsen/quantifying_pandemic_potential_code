library(tidyverse)
library(DescTools)
library(khroma)
library(MASS)
library(fitdistrplus)
library(ggpubr)
library(scales)

## function to add linear interpolations between measured datapoints
interpolation <- function(row1, row2, data, interval){
  index_1 <- data[row1, ]
  index_2 <- data[row2, ]
  times <- seq(index_1$dpe, index_2$dpe, interval)
  preds <- seq(index_1$nw_titer, index_2$nw_titer, length.out=length(times))
  df <- data.frame(dpe = times,
                   nw_titer = preds)
  return (df)
}

## create function for F(t)
calculate_pr_contact_pos <- function(lambda_integral){
  prob <- 1 - exp(-lambda_integral)
  return (prob)
}

interpolation_interval <- 0.001
LOD <- 0.5

## MLE for s
MLE_H1N1 <- 0.111
MLE_H3N2 <- 0.044

## define "contact" as one hour of exposure
exposure.length <- 1/24

## number of contacts per day (15 = rounded average for 5-9yo/day from Mossong et al.)
num.contacts <- 15

## number of trials
its <- 1000

plot_colors <- color("muted")(2)

set.seed(25)

# H1N1 analysis -----------------------------------------------------------

H1N1_ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H1N1_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "dpi", "nw_titer", "donor_dose")

H1N1_RC_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer) %>%
  ## removing the uninfected ferret
  filter(Ferret_ID != 7824) %>%
  mutate(dpe = dpi - 1)
H1N1_recipient_names <- unique(H1N1_RC_ferrets$Ferret_ID)

## generate interpolated viral loads for each contact ferret

H1N1_transmission_probs <- vector("list", length(H1N1_recipient_names))
names(H1N1_transmission_probs) <- H1N1_recipient_names

## generate predictions for each ferret
for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    dplyr::select(c(dpe, nw_titer))
  df_1.3 <- interpolation(1, 2, ferret_data, interpolation_interval)
  df_3.5 <- interpolation(2, 3, ferret_data, interpolation_interval)
  df_5.7 <- interpolation(3, 4, ferret_data, interpolation_interval)
  df_7.9 <- interpolation(4, 5, ferret_data, interpolation_interval)
  df_9.11 <- interpolation(5, 6, ferret_data, interpolation_interval)
  ## if infection resolves by 11 dpe
  if (df_9.11[[length(df_9.11$dpe), "nw_titer"]] == LOD){
    combo <- rbind(ferret_data, df_1.3, df_3.5, df_5.7, df_7.9, df_9.11)
    combo <- combo %>%
      ## get rid of duplicate rows
      distinct() %>%
      ## specify if values are interpolated or measured
      mutate(type = if_else(dpe %in% c(1, 3, 5, 7, 9, 11), "measured", "predicted")) %>%
      arrange(dpe) %>%
      ## ensure that times are numeric for future integration
      mutate(dpe = as.numeric(dpe))
  } else { ## if infection hasn't yet resolved, add another interpolation for assumed negative test
    times <- seq(df_9.11[[length(df_9.11$dpe), "dpe"]], 13, interpolation_interval)
    preds <- seq(df_9.11[[length(df_9.11$dpe), "nw_titer"]], 0.5, length.out=length(times))
    df_11.13 <- data.frame(dpe = times,
                           nw_titer = preds)
    combo <- rbind(ferret_data, df_1.3, df_3.5, df_5.7, df_7.9, df_9.11, df_11.13)
    combo <- combo %>%
      ## get rid of duplicate rows
      distinct() %>%
      ## specify if values are interpolated or measured
      mutate(type = if_else(dpe %in% c(1, 3, 5, 7, 9, 11), "measured", "predicted")) %>%
      arrange(dpe) %>%
      ## ensure that times are numeric for future integration
      mutate(dpe = as.numeric(dpe))
  }
  H1N1_transmission_probs[[ferret]] <- combo
}

## calculate probability of onwards transmission
## could collapse this step into the next loop; as-is, code caluclates prob of transmission
## for all time steps, which is not really necessary, but takes ~40sec on my laptop

for (ferret in H1N1_recipient_names){
  data <- H1N1_transmission_probs[[ferret]]
  lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(data$nw_titer)){
    if (data[[k,"nw_titer"]] == LOD){
      ## if titer = LOD, we assume that the force of infection is 0
      lambda <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
      lambda <- AUC(x=c(0, exposure.length), y=rep(data[k,"nw_titer"]*MLE_H1N1, 2), method="trapezoid")
    }
    lambda_vals[k] <- lambda
  }
  data$constant_foi <- lambda_vals
  ## compute the probability of transmission given this force of infection
  data$prob_transmission <- calculate_pr_contact_pos(data$constant_foi)
  H1N1_transmission_probs[[ferret]] <- data
}

## record individual R0s estimates for all trials
H1N1.indv.Z <- matrix(data=NA, ncol=its, nrow=length(H1N1_recipient_names))
## record negative binomial mu and k for each trial
H1N1.negb.fits <- matrix(data=NA, ncol=its, nrow=2)
rownames(H1N1.negb.fits) <- c("k", "mu")

## record individual gen time estimates for all trials
H1N1.gen.times <- matrix(data=NA, ncol=its, nrow=num.contacts*10)

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  num.offspring <- c()
  gen.time <- c()
  for (ferret in H1N1_recipient_names){
    data <- H1N1_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match
    contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=3)
    ## this keeps duplicate times
    pr_transmission <- c()
    for (t in contact.times){
      pr_transmission <- append(pr_transmission, data[[which(near(data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    random.draws <- runif(num.contacts*10)
    ## count number of infected contacts and record
    infected.contacts <- sum(pr_transmission > random.draws)
    num.offspring <- append(num.offspring, infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (data[[1,"nw_titer"]] > LOD){
      time.initial <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      time.initial <- data[[which.max(data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- contact.times[which(pr_transmission > random.draws)] - time.initial
    gen.time <- append(gen.time, gen.interval)
  }
  ## if offspring are generated, fit a distribution and add to tracker
  if (sum(num.offspring > 0)) {
    ## add offspring distribution to tracker
    H1N1.indv.Z[,i] <- num.offspring
    ## add neg B params to tracker
    H1N1.negb.fits[,i] <- fitdist(c(num.offspring), "nbinom", method="mle")$estimate
    ## add gen times to tracker
    H1N1.gen.times[1:length(gen.time), i] <- gen.time
    ## if no offspring are generated, only track Z
  } else {
    H1N1.indv.Z[,i] <- 0
    H1N1.negb.fits[2,i] <- 0
  }
}

## plot offspring distribution for one simulation

panel_a <- ggplot(as.data.frame(H1N1.indv.Z[,4]), aes(x=H1N1.indv.Z[, 4])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[1]]) +
  labs(x="Number of secondary cases", y="Number") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_c <- ggplot(as.data.frame(H1N1.gen.times[,4]), aes(x=H1N1.gen.times[, 4])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[1]]) +
  labs(x="Generation time (days)", y="Number") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

# H3N2 analysis -----------------------------------------------------------

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "dpi", "nw_titer", "donor_dose")

H3N2_RC_ferrets <- H3N2_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer) %>%
  ## exclude all negative tests at 13 and 14 dpi
  filter(dpi < 14) %>%
  mutate(dpe = dpi - 1)

H3N2_recipient_names <- unique(H3N2_RC_ferrets$Ferret_ID)

## exclude all animals that never became infected 
kept_names <- c()
for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  if (max(ferret_data$nw_titer) > LOD){
    kept_names <- append(kept_names, ferret)
  } else {}
}

H3N2_RC_ferrets <- H3N2_RC_ferrets %>%
  filter(Ferret_ID %in% kept_names)
H3N2_recipient_names <- unique(H3N2_RC_ferrets$Ferret_ID)

## generate interpolated viral loads for each contact ferret

H3N2_transmission_probs <- vector("list", length(H3N2_recipient_names))
names(H3N2_transmission_probs) <- H3N2_recipient_names

## generate predictions for each ferret
for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    dplyr::select(c(dpe, nw_titer))
  df_1.3 <- interpolation(1, 2, ferret_data, interpolation_interval)
  df_3.5 <- interpolation(2, 3, ferret_data, interpolation_interval)
  df_5.7 <- interpolation(3, 4, ferret_data, interpolation_interval)
  df_7.9 <- interpolation(4, 5, ferret_data, interpolation_interval)
  df_9.11 <- interpolation(5, 6, ferret_data, interpolation_interval)
  ## if infection resolves by 11 dpe
  if (df_9.11[[length(df_9.11$dpe), "nw_titer"]] == LOD){
    combo <- rbind(ferret_data, df_1.3, df_3.5, df_5.7, df_7.9, df_9.11)
    combo <- combo %>%
      ## get rid of duplicate rows
      distinct() %>%
      ## specify if values are interpolated or measured
      mutate(type = if_else(dpe %in% c(1, 3, 5, 7, 9, 11), "measured", "predicted")) %>%
      arrange(dpe) %>%
      ## ensure that times are numeric for future integration
      mutate(dpe = as.numeric(dpe))
  } else { ## if infection hasn't yet resolved, add another interpolation for assumed negative test
    times <- seq(df_9.11[[length(df_9.11$dpe), "dpe"]], 13, interpolation_interval)
    preds <- seq(df_9.11[[length(df_9.11$dpe), "nw_titer"]], 0.5, length.out=length(times))
    df_11.13 <- data.frame(dpe = times,
                           nw_titer = preds)
    combo <- rbind(ferret_data, df_1.3, df_3.5, df_5.7, df_7.9, df_9.11, df_11.13)
    combo <- combo %>%
      ## get rid of duplicate rows
      distinct() %>%
      ## specify if values are interpolated or measured
      mutate(type = if_else(dpe %in% c(1, 3, 5, 7, 9, 11), "measured", "predicted")) %>%
      arrange(dpe) %>%
      ## ensure that times are numeric for future integration
      mutate(dpe = as.numeric(dpe))
  }
  H3N2_transmission_probs[[ferret]] <- combo
}

## calculate probability of onwards transmission
## could collapse this step into the next loop; as-is, code caluclates prob of transmission
## for all time steps, which is not really necessary, but takes ~40sec on my laptop

for (ferret in H3N2_recipient_names){
  data <- H3N2_transmission_probs[[ferret]]
  lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(data$nw_titer)){
    if (data[[k,"nw_titer"]] == LOD){
      ## if titer = LOD, we assume that the force of infection is 0
      lambda <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
      lambda <- AUC(x=c(0, exposure.length), y=rep(data[k,"nw_titer"]*MLE_H3N2, 2), method="trapezoid")
    }
    lambda_vals[k] <- lambda
  }
  data$constant_foi <- lambda_vals
  ## compute the probability of transmission given this force of infection
  data$prob_transmission <- calculate_pr_contact_pos(data$constant_foi)
  H3N2_transmission_probs[[ferret]] <- data
}

## record individual R0s estimates for all trials
H3N2.indv.Z <- matrix(data=NA, ncol=its, nrow=length(H3N2_recipient_names))
## record negative binomial mu and k for each trial
H3N2.negb.fits <- matrix(data=NA, ncol=its, nrow=2)
rownames(H3N2.negb.fits) <- c("k", "mu")

## record individual gen time estimates for all trials
H3N2.gen.times <- matrix(data=NA, ncol=its, nrow=num.contacts*10)

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  num.offspring <- c()
  gen.time <- c()
  for (ferret in H3N2_recipient_names){
    data <- H3N2_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match
    contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=3)
    ## this keeps duplicate times
    pr_transmission <- c()
    for (t in contact.times){
      pr_transmission <- append(pr_transmission, data[[which(near(data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    random.draws <- runif(num.contacts*10)
    ## count number of infected contacts and record
    infected.contacts <- sum(pr_transmission > random.draws)
    num.offspring <- append(num.offspring, infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (data[[1,"nw_titer"]] > LOD){
      time.initial <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      time.initial <- data[[which.max(data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- contact.times[which(pr_transmission > random.draws)] - time.initial
    gen.time <- append(gen.time, gen.interval)
  }
  ## if offspring are generated, fit a distribution and add to tracker
  if (sum(num.offspring > 0)) {
    ## add offspring distribution to tracker
    H3N2.indv.Z[,i] <- num.offspring
    ## add neg B params to tracker
    H3N2.negb.fits[,i] <- fitdist(c(num.offspring), "nbinom", method="mle")$estimate
    ## add gen times to tracker
    H3N2.gen.times[1:length(gen.time), i] <- gen.time
    ## if no offspring are generated, only track Z
  } else {
    H3N2.indv.Z[,i] <- 0
    H3N2.negb.fits[2,i] <- 0
  }
}

## plot offspring distribution for one simulation

panel_b <- ggplot(as.data.frame(H3N2.indv.Z[,5]), aes(x=H3N2.indv.Z[, 5])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[2]]) +
  labs(x="Number of secondary cases", y="Number") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_d <- ggplot(as.data.frame(H3N2.gen.times[,5]), aes(x=H3N2.gen.times[, 5])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[2]]) +
  labs(x="Generation time (days)", y="Number") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

# joint plots -------------------------------------------------------------

## plot mu distribution and test significance 

mu.vals <- data.frame(mu = c(c(H1N1.negb.fits[2,]), c(H3N2.negb.fits[2,])), 
                      Subtype = c(rep("H1N1", its), rep("H3N2", its)))
H1N1.mean.mu <- mean(H1N1.negb.fits[2,])
H3N2.mean.mu <- mean(H3N2.negb.fits[2,])

t.test(H1N1.negb.fits[2,], H3N2.negb.fits[2,], alternative = "two.sided")

panel_e <- ggplot(mu.vals, aes(x=mu, fill=Subtype, color=Subtype)) +
  geom_density(alpha=0.7) +
  ## add lines for mean mu
  geom_vline(xintercept=H1N1.mean.mu, color=plot_colors[[1]], linewidth=2, linetype=2) +
  geom_vline(xintercept=H3N2.mean.mu, color=plot_colors[[2]], linewidth=2, linetype=2) +
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  guides(fill = guide_legend(override.aes = list(alpha=1))) +
  xlim(0, 4) +
  labs(x=expression(paste("Basic reproductive number ", R[0])), y="Density") +
  theme_light()

## plot k density curve

k.density <- data.frame(k = c(H1N1.negb.fits[1,], H3N2.negb.fits[1,]), 
                             Subtype = c(rep("H1N1", its), rep("H3N2", its)))

panel_f <- ggplot(k.density, aes(k, color=Subtype)) +
  stat_ecdf(geom="line", linewidth=2) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(limits=c(0, 2), breaks=c(0, 1, 2)) +
  labs(x="Overdispersion parameter k", y="Cumulative density") +
  theme_light()

## plot generation intervals and test significance

H1N1.Tc <- c(H1N1.gen.times)[!is.na(c(H1N1.gen.times))]
H3N2.Tc <- c(H3N2.gen.times)[!is.na(c(H3N2.gen.times))]
Tc.values <- data.frame(Tc = c(H1N1.Tc, H3N2.Tc), 
                        Subtype = c(rep("H1N1", length(H1N1.Tc)), rep("H3N2", length(H3N2.Tc))))

t.test(H1N1.Tc, H3N2.Tc, altnervative="two.sided")

panel_g <- ggplot(Tc.values, aes(x=Tc, color=Subtype)) +
  geom_density(linewidth=2) +
  ## add lines for mean gen time
  geom_vline(xintercept=mean(H1N1.Tc), color=plot_colors[[1]], linewidth=2, linetype=2) +
  geom_vline(xintercept=mean(H3N2.Tc), color=plot_colors[[2]], linewidth=2, linetype=2) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = seq(0, 12, 2), limits=c(0, 12)) +
  labs(x=expression(paste("Mean generation time ", T[c], " (days)")), y="Density") + 
  theme_light()

## all plots

top <- ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, labels=c("A", "B", "C", "D"), align="h")

bottom <- ggarrange(panel_e, panel_f, panel_g, ncol=3, labels=c("E", "F", "G"), align="h", common.legend = T, legend="right", widths = c(2, 1, 1))

ggarrange(top, bottom, nrow=2, align="v")
