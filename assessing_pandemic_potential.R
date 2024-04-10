library(tidyverse)
library(DescTools)
library(khroma)
library(MASS)

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

interpolation_interval <- 0.1
LOD <- 0.5

## MLE for s
MLE_H1N1 <- 0.111
MLE_H3N2 <- 0.044

## define "contact" as one hour of exposure
exposure.length <- 1/24

## number of contacts over their infection (rounded average for 5-9yo from Mossong et al.)
num.contacts <- 15

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

for (ferret in H1N1_recipient_names){
  data <- H1N1_transmission_probs[[ferret]]
  lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(data$nw_titer)){
    if (data[[k,"nw_titer"]] == LOD){
    ## if titer = LOD, we assume that the force of infection is 0
    lambda <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
    lambda <- AUC(x=c(0, exposure.length), y=rep(data[k,"nw_titer"], 2), method="trapezoid")
    }
    lambda_vals[k] <- lambda
  }
  data$constant_foi <- lambda_vals
  ## compute the probability of transmission given this force of infection
  data$prob_transmission <- calculate_pr_contact_pos(data$constant_foi)
  H1N1_transmission_probs[[ferret]] <- data
}

## record number of infected contacts for each ferret
H1N1.num.offspring <- c()
## record generation times for each infected contact for each ferret
H1N1.gen.time <- c()

for (ferret in H1N1_recipient_names){
  data <- H1N1_transmission_probs[[ferret]]
  ## draw random times for the contacts to occur
  contact.times <- round(x=runif(num.contacts, min=data[[1,"dpe"]], max=data[[length(data$dpe), "dpe"]]), digits=1)
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
  ## sanity check
  if (length(gen.time) != sum(num.offspring)){
    stop(paste("error at ID", ferret))
  }
}

## fit the simulated offspring frequencies to a negative binomial distribution 
## size = shape (k) and mean = mu
hist(num.offspring, right=F)
H1N1_negb_fit <- fitdistr(num.offspring, "negative binomial")
H1N1_mean_R0 <- H1N1_negb_fit$estimate[["mu"]]
H1N1_k <- H1N1_negb_fit$estimate[["size"]]

## plot probability distribution of fitted neg binom over histogram
x <- seq(0, num.contacts, 1)
negb_fit <- data.frame(x = x, 
                       y = dnbinom(x, size=H1N1_k, mu=H1N1_mean_R0))
ggplot(as.data.frame(num.offspring), aes(x=num.offspring)) +
  geom_histogram() +
  geom_line(data=negb_fit, aes(x=x, y=y)) +
  theme_light()
## plot just probability distribution as a barchart
ggplot(negb_fit, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  labs(x="Number of offspring", y="Probability") +
  theme_light()

## generation interval plot
ggplot(as.data.frame(gen.time), aes(x=gen.time)) +
  geom_histogram() +
  labs(x="Generation time (days)", y="Count") +
  theme_light()
## fit simulated generation intervals with gamma distribution
H1N1_gen_time_fit <- fitdistr(gen.time, "gamma")
H1N1_shape <- H1N1_gen_time_fit$estimate[["shape"]]
H1N1_rate <- H1N1_gen_time_fit$estimate[["rate"]]
## plot prob distribution of fitted gamma over histogram
gamma_fit <- data.frame(x = seq(0, 12, 0.1), 
                        y = dgamma(seq(0, 12, 0.1), shape=H1N1_shape, rate=H1N1_rate))
ggplot(as.data.frame(gen.time), aes(x=gen.time)) +
  geom_histogram() +
  geom_line(data=gamma_fit, aes(x=x, y=y)) +
  labs(x="Generation time (days)", y="Count") +
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

for (ferret in H3N2_recipient_names){
  data <- H3N2_transmission_probs[[ferret]]
  lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(data$nw_titer)){
    if (data[[k,"nw_titer"]] == LOD){
      ## if titer = LOD, we assume that the force of infection is 0
      lambda <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
      lambda <- AUC(x=c(0, exposure.length), y=rep(data[k,"nw_titer"], 2), method="trapezoid")
    }
    lambda_vals[k] <- lambda
  }
  data$constant_foi <- lambda_vals
  ## compute the probability of transmission given this force of infection
  data$prob_transmission <- calculate_pr_contact_pos(data$constant_foi)
  H3N2_transmission_probs[[ferret]] <- data
}

## record number of infected contacts for each ferret
H3N2.num.offspring <- c()
## record generation times for each infected contact for each ferret
H3N2.gen.time <- c()

for (ferret in H3N2_recipient_names){
  data <- H3N2_transmission_probs[[ferret]]
  ## draw random times for the contacts to occur
  contact.times <- round(x=runif(num.contacts, min=data[[1,"dpe"]], max=data[[length(data$dpe), "dpe"]]), digits=1)
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
  ## sanity check
  if (length(gen.time) != sum(num.offspring)){
    stop(paste("error at ID", ferret))
  }
}

## fit the simulated offspring frequencies to a negative binomial distribution 
## size = shape (k) and mean = mu
hist(num.offspring, right=F)
H3N2_negb_fit <- fitdistr(num.offspring, "negative binomial")
H3N2_mean_R0 <- H3N2_negb_fit$estimate[["mu"]]
H3N2_k <- H3N2_negb_fit$estimate[["size"]]

## plot probability distribution of fitted neg binom over histogram
x <- seq(0, num.contacts, 1)
negb_fit <- data.frame(x = x, 
                       y = dnbinom(x, size=H3N2_k, mu=H3N2_mean_R0))
ggplot(as.data.frame(num.offspring), aes(x=num.offspring)) +
  geom_histogram() +
  geom_line(data=negb_fit, aes(x=x, y=y)) +
  theme_light()
## plot just probability distribution as a barchart
ggplot(negb_fit, aes(x=x, y=y)) +
  geom_bar(stat="identity") +
  labs(x="Number of offspring", y="Probability") +
  theme_light()

## generation interval plot
ggplot(as.data.frame(gen.time), aes(x=gen.time)) +
  geom_histogram() +
  labs(x="Generation time (days)", y="Count") +
  theme_light()
## fit simulated generation intervals with gamma distribution
H3N2_gen_time_fit <- fitdistr(gen.time, "gamma")
H3N2_shape <- H3N2_gen_time_fit$estimate[["shape"]]
H3N2_rate <- H3N2_gen_time_fit$estimate[["rate"]]
## plot prob distribution of fitted gamma over histogram
gamma_fit <- data.frame(x = seq(0, 12, 0.1), 
                        y = dgamma(seq(0, 12, 0.1), shape=H3N2_shape, rate=H3N2_rate))
ggplot(as.data.frame(gen.time), aes(x=gen.time)) +
  geom_histogram() +
  geom_line(data=gamma_fit, aes(x=x, y=y)) +
  labs(x="Generation time (days)", y="Count") +
  theme_light()
