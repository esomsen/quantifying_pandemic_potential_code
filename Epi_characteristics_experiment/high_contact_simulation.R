library(tidyverse)
library(DescTools)
library(MASS)
library(fitdistrplus)

## function to add linear interpolations between measured datapoints
interpolation <- function(row1, row2, data, interval=0.01){
  tmp.times <- seq(data[[row1,"dpe"]], data[[row2,"dpe"]], interval)
  tmp.preds <- seq(data[[row1,"nw_titer"]], data[[row2,"nw_titer"]], length.out=length(tmp.times))
  tmp.df <- data.frame(dpe = tmp.times, nw_titer = tmp.preds)
  return (tmp.df)
  rm(list=ls(pattern="^tmp"))
}

## create function for F(t)
calculate_pr_contact_pos <- function(lambda_integral){
  prob <- 1 - exp(-lambda_integral)
  return (prob)
}

interpolation_interval <- 0.01
LOD <- 1

## MLE for s
MLE_H1N1 <- 0.111
MLE_H3N2 <- 0.047

## define "contact" as one hour of exposure
exposure.length <- 1/24

## number of contacts per day
num.contacts <- 1000

## number of trials
its <- 1000

plot_colors <- color("muted")(2)

set.seed(25)

# H1N1 analysis -----------------------------------------------------------

H1N1_ferrets <- read_csv("/home/esomsen/within-host/H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
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

for (ferret in H1N1_recipient_names){
  ## extract data for each donor
  tmp.data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    dplyr::select(c(dpe, nw_titer))
  ## create temporary df to store interpolations
  tmp.df <- data.frame()
  ## loop through all data until last timepoint
  for (t in 1:(length(tmp.data$dpe)-1)){
    ## interpolation
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data))
  }
  ## if no ending negative test, assume animal would have tested negative at next time
  if (tmp.df[length(tmp.df$dpe),"nw_titer"] > LOD){
    tmp.times <- seq(tmp.df[length(tmp.df$dpe),"dpe"], tmp.df[length(tmp.df$dpe),"dpe"]+2, interpolation_interval)
    tmp.preds <- seq(tmp.df[length(tmp.df$dpe),"nw_titer"], LOD, length.out=length(tmp.times))
    tmp.df <- rbind(tmp.df, data.frame(dpe = tmp.times, 
                                       nw_titer = tmp.preds))
  }
  ## force numeric, remove duplicate rows
  tmp.df <- tmp.df %>%
    distinct() %>%
    mutate(dpe = as.numeric(dpe))
  ## store interpolated titers in list
  H1N1_transmission_probs[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## calculate probability of onwards transmission
## could collapse this step into the next loop; as-is, code caluclates prob of transmission
## for all time steps, which is not really necessary, but very fast on my laptop

for (ferret in H1N1_recipient_names){
  tmp.data <- H1N1_transmission_probs[[ferret]]
  tmp.lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(tmp.data$nw_titer)){
    if (tmp.data[[k,"nw_titer"]] <= LOD){
      ## if titer <= LOD, we assume that the force of infection is 0
      tmp.lambda_vals[k] <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
      tmp.lambda_vals[k] <- AUC(x=c(0, exposure.length), y=rep(tmp.data[k,"nw_titer"]*MLE_H1N1, 2), method="trapezoid")
    }
  }
  tmp.data$constant_foi <- tmp.lambda_vals
  ## compute the probability of transmission given this force of infection
  tmp.data$prob_transmission <- calculate_pr_contact_pos(tmp.data$constant_foi)
  H1N1_transmission_probs[[ferret]] <- tmp.data
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## record individual R0s estimates for all trials
H1N1.indv.Z <- matrix(data=NA, ncol=its, nrow=length(H1N1_recipient_names))
## record negative binomial mu and k for each trial
H1N1.negb.fits <- matrix(data=NA, ncol=its, nrow=2)
rownames(H1N1.negb.fits) <- c("k", "mu")

## record individual gen time estimates for all trials
## 10 day experiment
H1N1.gen.times <- matrix(data=NA, ncol=its, nrow=num.contacts*10)

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  tmp.num.offspring <- c()
  tmp.gen.time <- c()
  for (ferret in H1N1_recipient_names){
    tmp.data <- H1N1_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## this keeps duplicate times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.data[[which(near(tmp.data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    tmp.draws <- runif(num.contacts*10)
    ## count number of infected contacts and record
    tmp.infected.contacts <- sum(tmp.pr.transmission > tmp.draws)
    tmp.num.offspring <- append(tmp.num.offspring, tmp.infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (tmp.data[[1,"nw_titer"]] > LOD){
      tmp.init.time <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      tmp.init.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- tmp.contact.times[which(tmp.pr.transmission > tmp.draws)] - tmp.init.time
    tmp.gen.time <- append(tmp.gen.time, gen.interval)
  }
  ## if offspring are generated, fit a distribution and add to tracker
  if (sum(tmp.num.offspring > 0)) {
    ## add offspring distribution to tracker
    H1N1.indv.Z[,i] <- tmp.num.offspring
    ## add neg B params to tracker
    H1N1.negb.fits[,i] <- fitdist(c(tmp.num.offspring), "nbinom", method="mle")$estimate
    ## add gen times to tracker
    H1N1.gen.times[1:length(tmp.gen.time), i] <- tmp.gen.time
    ## if no offspring are generated, only track Z
  } else {
    H1N1.indv.Z[,i] <- 0
    H1N1.negb.fits[2,i] <- 0
  }
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

# H3N2 analysis -----------------------------------------------------------

H3N2_ferrets <- read_csv("/home/esomsen/within-host/H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
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

for (ferret in H3N2_recipient_names){
  ## extract data for each donor
  tmp.data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    dplyr::select(c(dpe, nw_titer))
  ## create temporary df to store interpolations
  tmp.df <- data.frame()
  ## loop through all data until last timepoint
  for (t in 1:(length(tmp.data$dpe)-1)){
    ## interpolation
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data))
  }
  ## if no ending negative test, assume animal would have tested negative at next time
  if (tmp.df[length(tmp.df$dpe),"nw_titer"] > LOD){
    tmp.times <- seq(tmp.df[length(tmp.df$dpe),"dpe"], tmp.df[length(tmp.df$dpe),"dpe"]+2, interpolation_interval)
    tmp.preds <- seq(tmp.df[length(tmp.df$dpe),"nw_titer"], LOD, length.out=length(tmp.times))
    tmp.df <- rbind(tmp.df, data.frame(dpe = tmp.times, 
                                       nw_titer = tmp.preds))
  }
  ## force numeric, remove duplicate rows
  tmp.df <- tmp.df %>%
    distinct() %>%
    mutate(dpe = as.numeric(dpe))
  ## store interpolated titers in list
  H3N2_transmission_probs[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## calculate probability of onwards transmission
## could collapse this step into the next loop; as-is, code caluclates prob of transmission
## for all time steps, which is not really necessary, but very fast on my laptop

for (ferret in H3N2_recipient_names){
  tmp.data <- H3N2_transmission_probs[[ferret]]
  tmp.lambda_vals <- c()
  ## calculate a constant force of infection given 1 hour exposure to each viral titer
  for (k in 1:length(tmp.data$nw_titer)){
    if (tmp.data[[k,"nw_titer"]] <= LOD){
      ## if titer <= LOD, we assume that the force of infection is 0
      tmp.lambda_vals[k] <- 0 
    } else { ## otherwise, calculate constant lambda for one-hour exposure
      tmp.lambda_vals[k] <- AUC(x=c(0, exposure.length), y=rep(tmp.data[k,"nw_titer"]*MLE_H3N2, 2), method="trapezoid")
    }
  }
  tmp.data$constant_foi <- tmp.lambda_vals
  ## compute the probability of transmission given this force of infection
  tmp.data$prob_transmission <- calculate_pr_contact_pos(tmp.data$constant_foi)
  H3N2_transmission_probs[[ferret]] <- tmp.data
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## record individual R0s estimates for all trials
H3N2.indv.Z <- matrix(data=NA, ncol=its, nrow=length(H3N2_recipient_names))
## record negative binomial mu and k for each trial
H3N2.negb.fits <- matrix(data=NA, ncol=its, nrow=2)
rownames(H3N2.negb.fits) <- c("k", "mu")

## record individual gen time estimates for all trials
## 10 day experiment
H3N2.gen.times <- matrix(data=NA, ncol=its, nrow=num.contacts*10)

for (i in 1:its){
  ## record number of infected contacts for each ferret in each trial
  tmp.num.offspring <- c()
  tmp.gen.time <- c()
  for (ferret in H3N2_recipient_names){
    tmp.data <- H3N2_transmission_probs[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## this keeps duplicate times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.data[[which(near(tmp.data$dpe, t)),"prob_transmission"]])
    }
    ## draw a random number between 0-1; if transmission prob is higher than this number, contact is infected
    tmp.draws <- runif(num.contacts*10)
    ## count number of infected contacts and record
    tmp.infected.contacts <- sum(tmp.pr.transmission > tmp.draws)
    tmp.num.offspring <- append(tmp.num.offspring, tmp.infected.contacts)
    ## the generation interval is the time of successful infection - time infection begins in donor
    ## if the first test is > LOD, assume that the infection begins at 0dpe
    if (tmp.data[[1,"nw_titer"]] > LOD){
      tmp.init.time <- 0
    } else { ## if first test is negative, find the last time titer is at LOD
      tmp.init.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1,"dpe"]] 
    }
    ## for all timepoints at which contacts were infected, calculate the generation interval
    gen.interval <- tmp.contact.times[which(tmp.pr.transmission > tmp.draws)] - tmp.init.time
    tmp.gen.time <- append(tmp.gen.time, gen.interval)
  }
  ## if offspring are generated, fit a distribution and add to tracker
  if (sum(tmp.num.offspring > 0)) {
    ## add offspring distribution to tracker
    H3N2.indv.Z[,i] <- tmp.num.offspring
    ## add neg B params to tracker
    H3N2.negb.fits[,i] <- fitdist(c(tmp.num.offspring), "nbinom", method="mle")$estimate
    ## add gen times to tracker
    H3N2.gen.times[1:length(tmp.gen.time), i] <- tmp.gen.time
    ## if no offspring are generated, only track Z
  } else {
    H3N2.indv.Z[,i] <- 0
    H3N2.negb.fits[2,i] <- 0
  }
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

mean(c(H1N1.gen.times), na.rm=T)

mean(c(H3N2.gen.times), na.rm=T)

save.image(file="high_contact_sim.RData")