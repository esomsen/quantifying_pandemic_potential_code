library(tidyverse)
library(DescTools)
library(khroma)
library(MASS)
library(fitdistrplus)
library(ggpubr)

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

## number of contacts (rounded average for 5-9yo/day from Mossong et al.)
num.contacts <- 15

## number of trials
its <- 100

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

## record R0 and k values for all trials
H1N1.R0s <- c()
H1N1.ks <- c()
## record gamma shape and rate for all trials
H1N1.shapes <- c()
H1N1.rates <- c()
## record individual R0s estimates for all trials
H1N1.indv.Z <- matrix(data=0, ncol=its, nrow=length(H1N1_recipient_names))

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  H1N1.num.offspring <- c()
  ## record generation times for each infected contact for each ferret in each trial
  H1N1.gen.time <- c()
  for (ferret in H1N1_recipient_names){
    data <- H1N1_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur between 1-11dpe
    contact.times <- round(x=runif(num.contacts, min=1, max=11), digits=3)
    ## this keeps duplicate times
    pr_transmission <- c()
    for (t in contact.times){
      pr_transmission <- append(pr_transmission, data[[which(near(data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    random.draws <- runif(num.contacts)
    ## count number of infected contacts and record
    infected.contacts <- sum(pr_transmission > random.draws)
    H1N1.num.offspring <- append(H1N1.num.offspring, infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (data[[1,"nw_titer"]] > LOD){
      time.initial <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      time.initial <- data[[which.max(data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- contact.times[which(pr_transmission > random.draws)] - time.initial
    H1N1.gen.time <- append(H1N1.gen.time, gen.interval)
  }
  ## fit the simulated offspring frequencies to a negative binomial distribution 
  ## size = shape (k) and mean = mu
  H1N1_negb_fit <- fitdist(H1N1.num.offspring, "nbinom", method="mle")
  H1N1.R0s[i] <- H1N1_negb_fit$estimate[["mu"]]
  H1N1.ks[i] <- H1N1_negb_fit$estimate[["size"]]
  H1N1.indv.Z[,i] <- H1N1.num.offspring
  ## if there is only one offspring, use "delta function"
  if (sum(H1N1.num.offspring) == 1){
    ## I understand how to incorporate if I was tracking the mean (the mean would just be the estimate)
    ## but I'm not sure it's possible to relate the delta function to the shape and rate 
  } else {
    ## fit simulated generation intervals with gamma distribution 
    H1N1_gen_time_fit <- fitdist(H1N1.gen.time, "gamma", method="mle")
    H1N1.shapes[i] <- H1N1_gen_time_fit$estimate[["shape"]]
    H1N1.rates[i] <- H1N1_gen_time_fit$estimate[["rate"]] 
  }
}

H1N1_mean_R0 <- mean(H1N1.R0s)
H1N1_mean_k <- mean(H1N1.ks) ## is this the best way? Note that values are not normally dist
plot(H1N1.R0s, H1N1.ks)

H1N1_mean_shape <- mean(H1N1.shapes)
H1N1_mean_rate <- mean(H1N1.rates)

## save values for future use
save(H1N1_mean_R0, H1N1_mean_k, H1N1_mean_shape, H1N1_mean_rate, file="H1N1_distr.RData")

H1N1.indv.Z <- as.data.frame(H1N1.indv.Z)
H1N1.indv.Z$Ferret_ID <- as.character(H1N1_recipient_names)
H1N1.indv.Z <- H1N1.indv.Z %>%
  pivot_longer(cols=1:its, names_to="sim", names_prefix="V", values_to="Z")

## plot botplot of all ferret offspring
ggplot(H1N1.indv.Z, aes(x=Ferret_ID, y=Z)) +
  geom_boxplot() +
  stat_summary(fun="mean", geom="point", size=2) +
  labs(x="Ferret ID", y="Number of secondary cases") +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 

## plot average offspring distribution
H1N1_negb <- data.frame(x=seq(0, 10, 1), 
                        y=dnbinom(seq(0, 10, 1), size=H1N1_mean_k, mu=H1N1_mean_R0))
ggplot(H1N1_negb, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  labs(x="Number of offspring", y="Probability") +
  theme_light()

## plot generation interval
H1N1_gamma <- data.frame(x=seq(0, 12, 0.01), 
                         y=dgamma(seq(0, 12, 0.01), shape=H1N1_mean_shape, rate=H1N1_mean_rate))
ggplot(H1N1_gamma, aes(x=x, y=y)) +
  geom_line() +
  labs(x="Generation interval (days)", y="Probability") +
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

## record R0 and k values for all trials
H3N2.R0s <- c()
H3N2.ks <- c()
## record gamma shape and rate for all trials
H3N2.shapes <- c()
H3N2.rates <- c()
## record individual R0s estimates for all trials
H3N2.indv.Z <- matrix(data=0, ncol=its, nrow=length(H3N2_recipient_names))
#times <- c()

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  H3N2.num.offspring <- c()
  ## record generation times for each infected contact for each ferret in each trial
  H3N2.gen.time <- c()
  for (ferret in H3N2_recipient_names){
    data <- H3N2_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur between 1-11dpe
    contact.times <- round(x=runif(num.contacts, min=1, max=11), digits=3)
    ## this keeps duplicate times
    pr_transmission <- c()
    for (t in contact.times){
      pr_transmission <- append(pr_transmission, data[[which(near(data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    random.draws <- runif(num.contacts)
    ## count number of infected contacts and record
    infected.contacts <- sum(pr_transmission > random.draws)
    H3N2.num.offspring <- append(H3N2.num.offspring, infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (data[[1,"nw_titer"]] > LOD){
      time.initial <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      time.initial <- data[[which.max(data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- contact.times[which(pr_transmission > random.draws)] - time.initial
    H3N2.gen.time <- append(H3N2.gen.time, gen.interval)
  }
  ## fit the simulated offspring frequencies to a negative binomial distribution 
  ## size = shape (k) and mean = mu
  H3N2_negb_fit <- fitdist(H3N2.num.offspring, "nbinom", method="mle")
  H3N2.R0s[i] <- H3N2_negb_fit$estimate[["mu"]]
  H3N2.ks[i] <- H3N2_negb_fit$estimate[["size"]]
  H3N2.indv.Z[,i] <- H3N2.num.offspring
  ## fit simulated generation intervals with gamma distribution 
  #times <- append(times, H3N2.gen.time)
  H3N2_gen_time_fit <- fitdist(H3N2.gen.time, "gamma", method="mle")
  H3N2.shapes[i] <- H3N2_gen_time_fit$estimate[["shape"]]
  H3N2.rates[i] <- H3N2_gen_time_fit$estimate[["rate"]]
}

H3N2_mean_R0 <- mean(H3N2.R0s)
H3N2_mean_k <- mean(H3N2.ks) ## is this the best way? Note that values are not normally dist
plot(H3N2.R0s, H3N2.ks)

H3N2_mean_shape <- mean(H3N2.shapes)
H3N2_mean_rate <- mean(H3N2.rates)

## save values for future use
save(H3N2_mean_R0, H3N2_mean_k, H3N2_mean_shape, H3N2_mean_rate, file="H3N2_distr.RData")

H3N2.indv.Z <- as.data.frame(H3N2.indv.Z)
H3N2.indv.Z$Ferret_ID <- as.character(H3N2_recipient_names)
H3N2.indv.Z <- H3N2.indv.Z %>%
  pivot_longer(cols=1:its, names_to="sim", names_prefix="V", values_to="Z")

## plot botplot of all ferret offspring
ggplot(H3N2.indv.Z, aes(x=Ferret_ID, y=Z)) +
  geom_boxplot() +
  stat_summary(fun="mean", geom="point", size=2) +
  labs(x="Ferret ID", y="Number of secondary cases") +
  theme_light() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 

## plot average offspring distribution
H3N2_negb <- data.frame(x=seq(0, 10, 1), 
                        y=dnbinom(seq(0, 10, 1), size=H3N2_mean_k, mu=H3N2_mean_R0))
ggplot(H3N2_negb, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  labs(x="Number of offspring", y="Probability") +
  theme_light()

## plot generation interval
H3N2_gamma <- data.frame(x=seq(0, 12, 0.01), 
                         y=dgamma(seq(0, 12, 0.01), shape=H3N2_mean_shape, rate=H3N2_mean_rate))
ggplot(H3N2_gamma, aes(x=x, y=y)) +
  geom_line() +
  labs(x="Generation interval (days)", y="Probability") +
  theme_light()

# combined plots -------------------------------------------------------------------

plot_colors <- color("muted")(2)

## individual offspring numbers
indv.Z <- rbind(H1N1.indv.Z, H3N2.indv.Z)
indv.Z$virus <- c(rep("H1N1", length(H1N1.indv.Z$Z)), rep("H3N2", length(H3N2.indv.Z$Z)))

offspring.plot <- ggplot(indv.Z, aes(x=Ferret_ID, y=Z, color=virus)) +
  geom_boxplot() +
  stat_summary(fun="mean", geom="point", size=2) +
  labs(x="Ferret ID", y="Number of secondary cases") +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  #scale_x_discrete(label=c(H1N1_recipient_names, H3N2_recipient_names)) + 
  ## above does not reoder, just changes labels
  theme_light() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 

## offspring distribution
negb.fit <- data.frame(x=seq(0, 9, 1), 
                       H1N1=dnbinom(seq(0, 9, 1), size=H1N1_mean_k, mu=H1N1_mean_R0),
                       H3N2=dnbinom(seq(0, 9, 1), size=H3N2_mean_k, mu=H3N2_mean_R0))
negb.fit <- negb.fit %>%
              pivot_longer(cols=2:3, names_to="subtype", values_to="prob")

negb.plot <- ggplot(negb.fit, aes(x=x, y=prob, color=subtype, fill=subtype)) +
  geom_bar(stat = "identity", position = "dodge", alpha=0.7) +
  ## add lines for mean R0
  geom_vline(xintercept=H1N1_mean_R0, color=plot_colors[[1]], linewidth=2, linetype=2) +
  geom_vline(xintercept=H3N2_mean_R0, color=plot_colors[[2]], linewidth=2, linetype=2) +
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  guides(fill = guide_legend(override.aes = list(alpha=1))) +
  scale_x_continuous(breaks=seq(0, 8, 2)) +
  labs(x="Number of offspring", y="Probability") +
  theme_light()

## generation time
gamma.fit <- data.frame(x=seq(0, 12, 0.01),
                        H1N1=dgamma(seq(0, 12, 0.01), shape=H1N1_mean_shape, rate=H1N1_mean_rate),
                        H3N2=dgamma(seq(0, 12, 0.01), shape=H3N2_mean_shape, rate=H3N2_mean_rate))
gamma.fit <- gamma.fit %>%
  pivot_longer(cols=2:3, names_to="subtype", values_to="prob")

gen.plot <- ggplot(gamma.fit, aes(x=x, y=prob, color=subtype, fill=subtype)) +
  geom_line(linewidth=2) +
  #geom_vline(xintercept=H1N1_mean_gen_time, color=plot_colors[[1]], linewidth=2, linetype=2) +
  #geom_vline(xintercept=H3N2_mean_gen_time, color=plot_colors[[2]], linewidth=2, linetype=2) +
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]])) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_x_continuous(breaks=seq(0, 12, 2)) +
  labs(x="Generation time (days)", y="Probability") +
  theme_light()

ggarrange(offspring.plot, negb.plot, gen.plot, ncol=1, nrow=3)
