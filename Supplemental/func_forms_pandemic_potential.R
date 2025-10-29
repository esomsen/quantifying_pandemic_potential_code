library(tidyverse)
library(DescTools)
library(khroma)
library(MASS)
library(fitdistrplus)
library(ggpubr)

## function to add linear interpolations between measured datapoints
interpolation <- function(row1, row2, data, interval){
  fxn.times <- seq(data[[row1,"dpe"]], data[[row2,"dpe"]], interval)
  fxn.preds <- seq(data[[row1,"nw_titer"]], data[[row2,"nw_titer"]], length.out=length(fxn.times))
  fxn.df <- data.frame(dpe = fxn.times, nw_titer = fxn.preds)
  return (fxn.df)
  rm(list=ls(pattern="^fxn"))
}

## create function for F(t)
calculate_pr_contact_pos <- function(lambda_integral){
  prob <- 1 - exp(-lambda_integral)
  return (prob)
}

interpolation_interval <- 0.01
LOD <- 1

## define "contact" as one hour of exposure
exposure.length <- 1/24

## number of contacts per day (15 = rounded average for 5-9yo/day from Mossong et al.)
num.contacts <- 15

## number of trials
its <- 1000

plot_colors <- color("muted")(2)

set.seed(25)

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
H1N1_interpolated_titers <- vector("list", length(H1N1_recipient_names))
names(H1N1_interpolated_titers) <- H1N1_recipient_names

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
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data, interpolation_interval))
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
  H1N1_interpolated_titers[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

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
H3N2_interpolated_titers <- vector("list", length(H3N2_recipient_names))
names(H3N2_interpolated_titers) <- H3N2_recipient_names

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
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data, interpolation_interval))
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
  H3N2_interpolated_titers[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

# LINEAR ------------------------------------------------------------------

## probability trace for H1N1 and H3N2 from MLE
H1N1.MLE.logL <- -46.48756
H3N2.MLE.logL <- -60.16884
load(file="linear.prob.trace.Rdata")

## draw s values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H1N1.alt.s.vals <- c()
while (length(H1N1.alt.s.vals) < its){
  tmp.draw <- sample(H1N1.prob.trace$s, 1)
  tmp.accept.ratio <- exp(H1N1.prob.trace[which(near(H1N1.prob.trace$s, tmp.draw)),2] - H1N1.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H1N1.alt.s.vals <- append(H1N1.alt.s.vals, tmp.draw)
  }
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
    tmp.data <- H1N1_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= LOD){
        tmp.pr <- calculate_pr_contact_pos(AUC(x=c(0, exposure.length), 
                                               y=rep(10^tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"]*H1N1.alt.s.vals[i], 2), method="trapezoid"))
      } else { ## if titer < LOD, we assume that the force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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

## plot offspring distribution for one simulation

panel_a <- ggplot(as.data.frame(H1N1.indv.Z[,11]), aes(x=H1N1.indv.Z[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[1]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_c <- ggplot(as.data.frame(H1N1.gen.times[,11]), aes(x=H1N1.gen.times[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[1]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H1N1.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H1N1.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  tmp.gamma.dist <- fitdist(H1N1.gen.times[,i][!is.na(H1N1.gen.times[,i])], "gamma", "mle")
  ## mean
  H1N1.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
  ## coefficient of variation
  H1N1.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]])
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## draw s values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H3N2.alt.s.vals <- c()
while (length(H3N2.alt.s.vals) < its){
  tmp.draw <- sample(H3N2.prob.trace$s, 1)
  tmp.accept.ratio <- exp(H3N2.prob.trace[which(near(H3N2.prob.trace$s, tmp.draw)),2] - H3N2.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H3N2.alt.s.vals <- append(H3N2.alt.s.vals, tmp.draw)
  }
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
    tmp.data <- H3N2_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= LOD){
        tmp.pr <- calculate_pr_contact_pos(AUC(x=c(0, exposure.length), 
                                               y=rep(10^tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"]*H3N2.alt.s.vals[i], 2), method="trapezoid"))
      } else { ## if titer < LOD, we assume that the force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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


## plot offspring distribution for one simulation

panel_b <- ggplot(as.data.frame(H3N2.indv.Z[,14]), aes(x=H3N2.indv.Z[,14])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[2]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_d <- ggplot(as.data.frame(H3N2.gen.times[,14]), aes(x=H3N2.gen.times[,14])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[2]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H3N2.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H3N2.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) > 1){
    if (length(unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])) == 1){
      H3N2.gamma.fits[1,i] <- unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])
      H3N2.gamma.fits[2,i] <- 0
    } else{
      tmp.gamma.dist <- fitdist(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])], "gamma", "mle")
      ## mean
      H3N2.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
      ## coefficient of variation
      H3N2.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]]) 
    }
  } else if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) == 1){
    H3N2.gamma.fits[1,i] <- H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]
    H3N2.gamma.fits[2,i] <- 0
  }
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## plot mu distribution

mu.vals <- data.frame(mu = c(c(H1N1.negb.fits[2,]), c(H3N2.negb.fits[2,])), 
                      Virus = c(rep("H1N1", its), rep("H3N2", its)))

panel_e <- ggplot(mu.vals, aes(x=mu, fill=Virus, color=Virus)) +
  geom_density(alpha=0.7) +
  scale_fill_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  scale_color_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  xlim(0, 4) +
  labs(x=expression(paste("Basic reproduction number ", R[0])), y="Density") +
  theme_light() +
  theme(legend.position= "inside", legend.position.inside = c(0.8, 0.8), 
        legend.key.size = unit(0.8, "cm"), legend.title = element_text(size=12), legend.text = element_text(size=10))

## plot k cumulative density

k.vals <- matrix(data=NA, nrow=length(seq(0, 1.5, 0.01)), ncol=3)
k.vals[,1] <- seq(0, 1.5, 0.01)
for (i in 1:length(k.vals[,1])) {
  k.vals[i,2] <- length(which(H1N1.negb.fits[1, ] < k.vals[i,1])) 
  k.vals[i,3] <- length(which(H3N2.negb.fits[1, ] < k.vals[i,1])) 
}

## change to proportion 
k.vals[,2] <- k.vals[,2] / (its - length(which(is.na(H1N1.negb.fits[1,]))))
k.vals[,3] <- k.vals[,3] / (its - length(which(is.na(H3N2.negb.fits[1,]))))

## manipulate for plotting
k.vals <- as.data.frame(k.vals)
names(k.vals) <- c("k", "H1N1", "H3N2")
k.vals <- k.vals %>%
  pivot_longer(cols=2:3, names_to="Virus", values_to="prop")

panel_f <- ggplot(k.vals, aes(x=k, y=prop, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = plot_colors) +
  labs(x="Overdispersion parameter k", y="Cumulative density") +
  ylim(0, 0.25) +
  theme_light() +
  theme(legend.position = "none")

## plot scatter of gamma mean and cov
gamma.mean.cov <- data.frame(mean = c(H1N1.gamma.fits[1,], H3N2.gamma.fits[1,]),
                             cov = c(H1N1.gamma.fits[2,], H3N2.gamma.fits[2,]), 
                             Virus = c(rep("H1N1", ncol(H1N1.gamma.fits)), rep("H3N2", ncol(H3N2.gamma.fits))))

panel_g <- ggplot(gamma.mean.cov, aes(x=mean, y=cov, color=Virus)) +
  geom_point(size=1, alpha=0.7) +
  geom_point(aes(x=2.85, y=0.30), color="#44AA99", size=2) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits=c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5), limits=c(0, 1.5)) +
  labs(x=expression(paste("Mean generation time ", T[c], " (days)")), y="Coefficient of variation") + 
  theme_light() +
  theme(legend.position = "none")

## all plots

top <- ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, labels=c("A", "B", "C", "D"), align="h")

bottom <- ggarrange(panel_e, panel_f, panel_g, ncol=3, labels=c("E", "F", "G"), align="h", common.legend = F, widths = c(2, 1, 1))

ggarrange(top, bottom, nrow=2, align="v")


# THRESHOLD -----------------------------------------------------------

## probability trace for H1N1 and H3N2 from MLE
H1N1.MLE.logL <- -16.87404
H3N2.MLE.logL <- -29.1636
load(file="threshold.prob.trace.Rdata")

## draw new param values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H1N1.alt.params <- matrix(data=NA, nrow=its, ncol=2)
colnames(H1N1.alt.params) <- c("h", "s")
while (is.na(H1N1.alt.params[its,1]) == TRUE){
  tmp.prop.h <- sample(as.numeric(rownames(H1N1.joint.log.probs)), 1)
  tmp.prop.s <- sample(as.numeric(colnames(H1N1.joint.log.probs)), 1)
  tmp.prop.logL <- H1N1.joint.log.probs[which(rownames(H1N1.joint.log.probs) == tmp.prop.h), which(colnames(H1N1.joint.log.probs) == tmp.prop.s)]
  tmp.accept.ratio <- exp(tmp.prop.logL - H1N1.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H1N1.alt.params[which.min(!is.na(H1N1.alt.params[,1])),1] <- tmp.prop.h
    H1N1.alt.params[which.min(!is.na(H1N1.alt.params[,1])),2] <- tmp.prop.s
  }
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
    tmp.data <- H1N1_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      ## if titer is above h, FOI = s (aka r in manuscript)
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= H1N1.alt.params[[i,1]]){
        tmp.pr <- calculate_pr_contact_pos(AUC(x=c(0, exposure.length), y=rep(H1N1.alt.params[[i,2]], 2), method="trapezoid"))
      } else { ## if titer < h, force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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

## plot offspring distribution for one simulation

panel_a <- ggplot(as.data.frame(H1N1.indv.Z[,11]), aes(x=H1N1.indv.Z[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[1]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_c <- ggplot(as.data.frame(H1N1.gen.times[,11]), aes(x=H1N1.gen.times[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[1]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H1N1.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H1N1.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  tmp.gamma.dist <- fitdist(H1N1.gen.times[,i][!is.na(H1N1.gen.times[,i])], "gamma", "mle")
  ## mean
  H1N1.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
  ## coefficient of variation
  H1N1.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]])
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## draw new param values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H3N2.alt.params <- matrix(data=NA, nrow=its, ncol=2)
colnames(H3N2.alt.params) <- c("h", "s")
while (is.na(H3N2.alt.params[its,1]) == TRUE){
  tmp.prop.h <- sample(as.numeric(rownames(H3N2.joint.log.probs)), 1)
  tmp.prop.s <- sample(as.numeric(colnames(H3N2.joint.log.probs)), 1)
  tmp.prop.logL <- H3N2.joint.log.probs[which(rownames(H3N2.joint.log.probs) == tmp.prop.h), which(colnames(H3N2.joint.log.probs) == tmp.prop.s)]
  tmp.accept.ratio <- exp(tmp.prop.logL - H3N2.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H3N2.alt.params[which.min(!is.na(H3N2.alt.params[,1])),1] <- tmp.prop.h
    H3N2.alt.params[which.min(!is.na(H3N2.alt.params[,1])),2] <- tmp.prop.s
  }
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
    tmp.data <- H3N2_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      ## if titer is above h, FOI = s (aka r in manuscript)
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= H3N2.alt.params[[i,1]]){
        tmp.pr <- calculate_pr_contact_pos(AUC(x=c(0, exposure.length), y=rep(H3N2.alt.params[[i,2]], 2), method="trapezoid"))
      } else { ## if titer < h, force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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

## plot offspring distribution for one simulation

panel_b <- ggplot(as.data.frame(H3N2.indv.Z[,4]), aes(x=H3N2.indv.Z[,4])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[2]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_d <- ggplot(as.data.frame(H3N2.gen.times[,4]), aes(x=H3N2.gen.times[,4])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[2]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H3N2.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H3N2.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) > 1){
    if (length(unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])) == 1){
      H3N2.gamma.fits[1,i] <- unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])
      H3N2.gamma.fits[2,i] <- 0
    } else{
      tmp.gamma.dist <- fitdist(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])], "gamma", "mle")
      ## mean
      H3N2.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
      ## coefficient of variation
      H3N2.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]]) 
    }
  } else if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) == 1){
    H3N2.gamma.fits[1,i] <- H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]
    H3N2.gamma.fits[2,i] <- 0
  }
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## plot mu distribution

mu.vals <- data.frame(mu = c(c(H1N1.negb.fits[2,]), c(H3N2.negb.fits[2,])), 
                      Virus = c(rep("H1N1", its), rep("H3N2", its)))

panel_e <- ggplot(mu.vals, aes(x=mu, fill=Virus, color=Virus)) +
  geom_density(alpha=0.7) +
  scale_fill_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  scale_color_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  xlim(0, 4) +
  labs(x=expression(paste("Basic reproduction number ", R[0])), y="Density") +
  theme_light() +
  theme(legend.position= "inside", legend.position.inside = c(0.8, 0.8), 
        legend.key.size = unit(0.8, "cm"), legend.title = element_text(size=12), legend.text = element_text(size=10))

## plot k cumulative density

k.vals <- matrix(data=NA, nrow=length(seq(0, 1.5, 0.01)), ncol=3)
k.vals[,1] <- seq(0, 1.5, 0.01)
for (i in 1:length(k.vals[,1])) {
  k.vals[i,2] <- length(which(H1N1.negb.fits[1, ] < k.vals[i,1])) 
  k.vals[i,3] <- length(which(H3N2.negb.fits[1, ] < k.vals[i,1])) 
}

## change to proportion 
k.vals[,2] <- k.vals[,2] / (its - length(which(is.na(H1N1.negb.fits[1,]))))
k.vals[,3] <- k.vals[,3] / (its - length(which(is.na(H3N2.negb.fits[1,]))))

## manipulate for plotting
k.vals <- as.data.frame(k.vals)
names(k.vals) <- c("k", "H1N1", "H3N2")
k.vals <- k.vals %>%
  pivot_longer(cols=2:3, names_to="Virus", values_to="prop")

panel_f <- ggplot(k.vals, aes(x=k, y=prop, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = plot_colors) +
  labs(x="Overdispersion parameter k", y="Cumulative density") +
  ylim(0, 0.25) +
  theme_light() +
  theme(legend.position = "none")

## plot scatter of gamma mean and cov
gamma.mean.cov <- data.frame(mean = c(H1N1.gamma.fits[1,], H3N2.gamma.fits[1,]),
                             cov = c(H1N1.gamma.fits[2,], H3N2.gamma.fits[2,]), 
                             Virus = c(rep("H1N1", ncol(H1N1.gamma.fits)), rep("H3N2", ncol(H3N2.gamma.fits))))

panel_g <- ggplot(gamma.mean.cov, aes(x=mean, y=cov, color=Virus)) +
  geom_point(size=1, alpha=0.7) +
  geom_point(aes(x=3.57, y=0.47), color="#44AA99", size=2) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits=c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5), limits=c(0, 1.5)) +
  labs(x=expression(paste("Mean generation time ", T[c], " (days)")), y="Coefficient of variation") + 
  theme_light() +
  theme(legend.position = "none")

## all plots

top <- ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, labels=c("A", "B", "C", "D"), align="h")

bottom <- ggarrange(panel_e, panel_f, panel_g, ncol=3, labels=c("E", "F", "G"), align="h", common.legend = F, widths = c(2, 1, 1))

ggarrange(top, bottom, nrow=2, align="v")


# HILL --------------------------------------------------------------------

## probability trace for H1N1 and H3N2 from MLE
H1N1.MLE.logL <- -18.93822
H3N2.MLE.logL <- -31.09032
load(file="H1.hill.prob.trace.RData")
load(file="H3.hill.prob.trace.RData")

## draw new param values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H1N1.alt.params <- matrix(data=NA, nrow=its, ncol=3)
colnames(H1N1.alt.params) <- c("q", "ka", "n")
while (is.na(H1N1.alt.params[its,1]) == TRUE){
  tmp.prop.q <- sample(as.numeric(dimnames(H1N1.log.probs)[[1]]), 1)
  tmp.prop.ka <- sample(as.numeric(dimnames(H1N1.log.probs)[[2]]), 1)
  tmp.prop.n <- sample(as.numeric(dimnames(H1N1.log.probs)[[3]]), 1)
  tmp.prop.logL <- H1N1.log.probs[which(dimnames(H1N1.log.probs)[[1]] == tmp.prop.q), which(dimnames(H1N1.log.probs)[[2]] == tmp.prop.ka), which(dimnames(H1N1.log.probs)[[3]] == tmp.prop.n)]
  tmp.accept.ratio <- exp(tmp.prop.logL - H1N1.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H1N1.alt.params[which.min(!is.na(H1N1.alt.params[,1])),1] <- tmp.prop.q
    H1N1.alt.params[which.min(!is.na(H1N1.alt.params[,1])),2] <- tmp.prop.ka
    H1N1.alt.params[which.min(!is.na(H1N1.alt.params[,1])),3] <- tmp.prop.n
  }
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
    tmp.data <- H1N1_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      ## if titer is above LOD, use Hill eq
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= LOD){
        ## Hill function
        tmp.theta <- (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"])^H1N1.alt.params[i,3] / (H1N1.alt.params[i,2]^H1N1.alt.params[i,3] + (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"])^H1N1.alt.params[i,3])
        tmp.pr <- AUC(x=c(0, exposure.length), rep(tmp.theta*H1N1.alt.params[i,1],2), method="trapezoid")
      } else { ## if titer < LOD, force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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


## plot offspring distribution for one simulation

panel_a <- ggplot(as.data.frame(H1N1.indv.Z[,11]), aes(x=H1N1.indv.Z[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[1]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_c <- ggplot(as.data.frame(H1N1.gen.times[,11]), aes(x=H1N1.gen.times[, 11])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[1]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H1N1.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H1N1.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  tmp.gamma.dist <- fitdist(H1N1.gen.times[,i][!is.na(H1N1.gen.times[,i])], "gamma", "mle")
  ## mean
  H1N1.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
  ## coefficient of variation
  H1N1.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]])
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

### H3N2

## draw new param values for each simulation
## use Metropolis-Hastings algorithm to accept or reject draws
H3N2.alt.params <- matrix(data=NA, nrow=its, ncol=3)
colnames(H3N2.alt.params) <- c("q", "ka", "n")
while (is.na(H3N2.alt.params[its,1]) == TRUE){
  tmp.prop.q <- sample(as.numeric(dimnames(H3N2.log.probs)[[1]]), 1)
  tmp.prop.ka <- sample(as.numeric(dimnames(H3N2.log.probs)[[2]]), 1)
  tmp.prop.n <- sample(as.numeric(dimnames(H3N2.log.probs)[[3]]), 1)
  tmp.prop.logL <- H3N2.log.probs[which(dimnames(H3N2.log.probs)[[1]] == tmp.prop.q), which(dimnames(H3N2.log.probs)[[2]] == tmp.prop.ka), which(dimnames(H3N2.log.probs)[[3]] == tmp.prop.n)]
  tmp.accept.ratio <- exp(tmp.prop.logL - H3N2.MLE.logL)
  if (tmp.accept.ratio > runif(1, 0, 1)){
    H3N2.alt.params[which.min(!is.na(H3N2.alt.params[,1])),1] <- tmp.prop.q
    H3N2.alt.params[which.min(!is.na(H3N2.alt.params[,1])),2] <- tmp.prop.ka
    H3N2.alt.params[which.min(!is.na(H3N2.alt.params[,1])),3] <- tmp.prop.n
  }
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
    tmp.data <- H3N2_interpolated_titers[[ferret]]
    ## draw random times for the contacts to occur each day, and round to match time resolution we have
    tmp.contact.times <- round(x=runif(num.contacts*10, min=seq(1, 10), max=seq(2, 11)), digits=2)
    ## calculate probability of transmission at these times
    tmp.pr.transmission <- c()
    for (t in tmp.contact.times){
      ## calculate a constant force of infection given 1 hour exposure
      ## if titer is above LOD, use Hill eq
      if (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"] >= LOD){
        ## Hill function
        tmp.theta <- (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"])^H3N2.alt.params[i,3] / (H3N2.alt.params[i,2]^H3N2.alt.params[i,3] + (tmp.data[which(near(tmp.data$dpe, t)),"nw_titer"])^H3N2.alt.params[i,3])
        tmp.pr <- AUC(x=c(0, exposure.length), rep(tmp.theta*H3N2.alt.params[i,1],2), method="trapezoid")
      } else { ## if titer < LOD, force of infection is 0
        tmp.pr <- 0
      }
      tmp.pr.transmission <- append(tmp.pr.transmission, tmp.pr)
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
  if (sum(tmp.num.offspring) > 0) {
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


## plot offspring distribution for one simulation

panel_b <- ggplot(as.data.frame(H3N2.indv.Z[,4]), aes(x=H3N2.indv.Z[,4])) +
  geom_histogram(closed="right", center=1, binwidth=0.5, fill=plot_colors[[2]]) +
  labs(x="Number of secondary cases", y="Count") +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5), limits=c(-0.5, 5.5)) +
  ylim(0, 7) +
  theme_light()

## plot generation interval distribution for one simulation

panel_d <- ggplot(as.data.frame(H3N2.gen.times[,4]), aes(x=H3N2.gen.times[,4])) +
  geom_histogram(closed="right", center=1, binwidth=0.25, fill=plot_colors[[2]]) +
  labs(x="Generation time (days)", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(-0.5, 10.5)) +
  ylim(0, 5) +
  theme_light()

## find number of simulations with k <= 1

length(which(H3N2.negb.fits[1,] <= 1))

## fit gamma distribution to all generation time estimates for each it
H3N2.gamma.fits <- matrix(data=NA, nrow=2, ncol=its)
for (i in 1:its){
  if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) > 1){
    if (length(unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])) == 1){
      H3N2.gamma.fits[1,i] <- unique(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])])
      H3N2.gamma.fits[2,i] <- 0
    } else{
      tmp.gamma.dist <- fitdist(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])], "gamma", "mle")
      ## mean
      H3N2.gamma.fits[1,i] <- tmp.gamma.dist$estimate[["shape"]] / tmp.gamma.dist$estimate[["rate"]]
      ## coefficient of variation
      H3N2.gamma.fits[2,i] <- 1 / sqrt(tmp.gamma.dist$estimate[["shape"]]) 
    }
  } else if (length(H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]) == 1){
    H3N2.gamma.fits[1,i] <- H3N2.gen.times[,i][!is.na(H3N2.gen.times[,i])]
    H3N2.gamma.fits[2,i] <- 0
  }
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## plot mu distribution 

mu.vals <- data.frame(mu = c(c(H1N1.negb.fits[2,]), c(H3N2.negb.fits[2,])), 
                      Virus = c(rep("H1N1", its), rep("H3N2", its)))

panel_e <- ggplot(mu.vals, aes(x=mu, fill=Virus, color=Virus)) +
  geom_density(alpha=0.7) +
  scale_fill_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  scale_color_manual(values = plot_colors, labels=c("Cal/2009", "Hong Kong/1968")) +
  xlim(0, 4) +
  labs(x=expression(paste("Basic reproduction number ", R[0])), y="Density") +
  theme_light() +
  theme(legend.position= "inside", legend.position.inside = c(0.8, 0.8), 
        legend.key.size = unit(0.8, "cm"), legend.title = element_text(size=12), legend.text = element_text(size=10))

## plot k cumulative density

k.vals <- matrix(data=NA, nrow=length(seq(0, 1.5, 0.01)), ncol=3)
k.vals[,1] <- seq(0, 1.5, 0.01)
for (i in 1:length(k.vals[,1])) {
  k.vals[i,2] <- length(which(H1N1.negb.fits[1, ] < k.vals[i,1])) 
  k.vals[i,3] <- length(which(H3N2.negb.fits[1, ] < k.vals[i,1])) 
}

## change to proportion 
k.vals[,2] <- k.vals[,2] / (its - length(which(is.na(H1N1.negb.fits[1,]))))
k.vals[,3] <- k.vals[,3] / (its - length(which(is.na(H3N2.negb.fits[1,]))))

## manipulate for plotting
k.vals <- as.data.frame(k.vals)
names(k.vals) <- c("k", "H1N1", "H3N2")
k.vals <- k.vals %>%
  pivot_longer(cols=2:3, names_to="Virus", values_to="prop")

panel_f <- ggplot(k.vals, aes(x=k, y=prop, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = plot_colors) +
  labs(x="Overdispersion parameter k", y="Cumulative density") +
  ylim(0, 0.25) +
  theme_light() +
  theme(legend.position = "none")

## plot scatter of gamma mean and cov
gamma.mean.cov <- data.frame(mean = c(H1N1.gamma.fits[1,], H3N2.gamma.fits[1,]),
                             cov = c(H1N1.gamma.fits[2,], H3N2.gamma.fits[2,]), 
                             Virus = c(rep("H1N1", ncol(H1N1.gamma.fits)), rep("H3N2", ncol(H3N2.gamma.fits))))

panel_g <- ggplot(gamma.mean.cov, aes(x=mean, y=cov, color=Virus)) +
  geom_point(size=1, alpha=0.7) +
  geom_point(aes(x=3.495, y=0.623), color="#44AA99", size=2) +
  scale_color_manual(values = plot_colors) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits=c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5), limits=c(0, 1.5)) +
  labs(x=expression(paste("Mean generation time ", T[c], " (days)")), y="Coefficient of variation") + 
  theme_light() +
  theme(legend.position = "none")

## all plots

top <- ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, labels=c("A", "B", "C", "D"), align="h")

bottom <- ggarrange(panel_e, panel_f, panel_g, ncol=3, labels=c("E", "F", "G"), align="h", common.legend = F, widths = c(2, 1, 1))

ggarrange(top, bottom, nrow=2, align="v")
