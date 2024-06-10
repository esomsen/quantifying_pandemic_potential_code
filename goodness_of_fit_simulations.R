library(tidyverse)
library(khroma)
library(DescTools)
library(ggpubr)

interpolation <- function(row1, row2, data, interval){
  index_1 <- data[row1, ]
  index_2 <- data[row2, ]
  times <- seq(index_1$dpe, index_2$dpe, interval)
  preds <- seq(index_1$nw_titer, index_2$nw_titer, length.out=length(times))
  df <- data.frame(dpe = times,
                   nw_titer = preds)
  return (df)
}

## resolution of interpolation
interpolation_interval <- 0.001
## bin length 
bin.length <- 1/24
## number of simulations
n.sims <- 1000

LOD <- 0.5

# H1N1 --------------------------------------------------------------------

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[[1]]

H1N1_ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H1N1_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

H1N1_DI_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "DI") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, dose, dpi, nw_titer) %>%
  mutate(dpe = dpi - 1)
H1N1_donor_names <- unique(H1N1_DI_ferrets$Ferret_ID)

H1N1_RC_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer, DI_RC_Pair) %>%
  ## removing the uninfected ferret
  filter(Ferret_ID != 7824) %>%
  mutate(dpe = dpi - 1)
H1N1_recipient_names <- unique(H1N1_RC_ferrets$Ferret_ID)

H1N1.MLE <- 0.111

H1N1.actual.or <- 37.703

## transmission simulations

H1N1_transmission_probs <- vector("list", length(H1N1_donor_names))
names(H1N1_transmission_probs) <- H1N1_donor_names

## generate interpolations for each ferret
for (ferret in H1N1_donor_names){
  ferret_data <- H1N1_DI_ferrets %>%
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

H1N1.odds.ratio <- c()

for (n in 1:n.sims){
  prob.transmission <- data.frame()
  ## find probability of transmission over entire infection
  for (ferret in H1N1_donor_names){
    ferret_data <- H1N1_transmission_probs[[ferret]]
    ## fix AUC vs lambda
    AUC <- AUC(x=ferret_data$dpe, y=ferret_data$nw_titer - LOD, method="trapezoid")
    prob.transmission <- rbind(prob.transmission, c(ferret, AUC, 1 - exp(-AUC*H1N1.MLE)))
  }
  names(prob.transmission) <- c("Ferret_ID", "AUC", "prob")
  ## stochastic draws
  prob.transmission$random.draws <- runif(length(H1N1_donor_names), 0, 1)
  prob.transmission$outcome <- ifelse(prob.transmission$prob > prob.transmission$random.draws, 1, 0)
  ## logistic regression for infection outcome by AUC
  logit <- glm(outcome ~ as.numeric(AUC), data=prob.transmission, family="binomial")
  ## the log odds of AUC means that for each unit increase in AUC, the odds a contact
  ## becomes infected increases by a factor = exp(coef(logit AUC))
  H1N1.odds.ratio[n] <- exp(logit$coefficients[[2]])
}

## for now, remove entries from when algorithm doesn't converge (perfect separation?)



## forward start time simulations

H1N1.simulated.start.times <- matrix(data=0, ncol=length(H1N1_donor_names), nrow=n.sims)
colnames(H1N1.simulated.start.times) <- H1N1_donor_names

for (ferret in H1N1_donor_names){
  data <- H1N1_transmission_probs[[ferret]]
  bins <- round(seq(from=data[[1,1]], to=data[[length(data$dpe),1]], by=bin.length), digits=3)
  lambdas <- data %>%
    filter(dpe %in% bins) %>%
    ## no transmission if titer is at LOD
    filter(nw_titer > LOD) %>%
    pull(nw_titer) * H1N1.MLE
  ## do I need to multiply by time?
  transmission.probs <- 1 - exp(-lambdas*bin.length)
  ## stochastic simulations
  for (n in 1:n.sims){
    random.draws <- runif(length(transmission.probs), 0, 1)
    ## what to do if there is no transmission?
    ## earliest simulated start time
    H1N1.simulated.start.times[n, ferret] <- bins[which(transmission.probs > random.draws)[1]]
  }
}

H1N1.simulated.start.times <- as.data.frame(H1N1.simulated.start.times)
H1N1.simulated.start.times$sim <- 1:n.sims
H1N1.simulated.start.times <- H1N1.simulated.start.times %>%
  pivot_longer(cols=1:length(H1N1_donor_names), names_to="animal", values_to="start.time")

## find first positive test for all contact ferrets
H1N1.first.positive.test <- data.frame()

for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  H1N1.first.positive.test <- rbind(H1N1.first.positive.test, ferret_data[which.max(ferret_data$nw_titer > LOD),c("DI_RC_Pair","dpe")])
}
names(H1N1.first.positive.test) <- c("animal", "start.time")
H1N1.first.positive.test$animal <- as.character(H1N1.first.positive.test$animal)

panel_b <- ggplot() +
  geom_violin(data=H1N1.simulated.start.times, aes(x=animal, y=start.time)) +
  geom_point(data=H1N1.first.positive.test, aes(x=animal, y=start.time), color=H1N1_color) +
  labs(x="Index", y="Infection start time") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2)) +
  theme_light()

# H3N2 --------------------------------------------------------------------

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[[2]]

LOD <- 0.5

H3N2.MLE <- 0.044

H3N2.actual.or <- 1.170589

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

H3N2_DI_ferrets <- H3N2_ferrets %>%
  filter(DI_RC == "DI") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, dose, dpi, nw_titer) %>%
  ## exclude all negative tests at 13 and 14 dpi
  filter(dpi < 14) %>%
  mutate(dpe = dpi - 1)

H3N2_donor_names <- unique(H3N2_DI_ferrets$Ferret_ID)

## exclude all animals that never became infected 
kept_names <- c()
for (ferret in H3N2_donor_names){
  ferret_data <- H3N2_DI_ferrets %>%
    filter(Ferret_ID == ferret)
  if (max(ferret_data$nw_titer, na.rm=T) > LOD){
    kept_names <- append(kept_names, ferret)
  } else {}
}

H3N2_DI_ferrets <- H3N2_DI_ferrets %>%
  filter(Ferret_ID %in% kept_names)
H3N2_donor_names <- unique(H3N2_DI_ferrets$Ferret_ID)

H3N2_RC_ferrets <- H3N2_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer, DI_RC_Pair) %>%
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

## transmission simulations

## exclude F6335 because of missing titer value

H3N2_donor_names <- H3N2_donor_names[c(1:7, 9:19)]
H3N2_DI_ferrets <- H3N2_DI_ferrets %>%
  filter(Ferret_ID %in% H3N2_donor_names)

H3N2_transmission_probs <- vector("list", length(H3N2_donor_names))
names(H3N2_transmission_probs) <- H3N2_donor_names

## generate interpolations for each ferret
for (ferret in H3N2_donor_names){
  ferret_data <- H3N2_DI_ferrets %>%
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

H3N2.odds.ratio <- c()

for (n in 1:n.sims){
  prob.transmission <- data.frame()
  ## find probability of transmission over entire infection
  for (ferret in H3N2_donor_names){
    ferret_data <- H3N2_transmission_probs[[ferret]]
    ## fix AUC vs lambda
    AUC <- AUC(x=ferret_data$dpe, y=ferret_data$nw_titer - LOD, method="trapezoid")
    prob.transmission <- rbind(prob.transmission, c(ferret, AUC, 1 - exp(-AUC*H3N2.MLE)))
  }
  names(prob.transmission) <- c("Ferret_ID", "AUC", "prob")
  ## stochastic draws
  prob.transmission$random.draws <- runif(length(H3N2_donor_names), 0, 1)
  prob.transmission$outcome <- ifelse(prob.transmission$prob > prob.transmission$random.draws, 1, 0)
  ## logistic regression for infection outcome by AUC
  logit <- glm(outcome ~ as.numeric(AUC), data=prob.transmission, family="binomial")
  ## the log odds of AUC means that for each unit increase in AUC, the odds a contact
  ## becomes infected increases by a factor = exp(coef(logit AUC))
  H3N2.odds.ratio[n] <- exp(logit$coefficients[[2]])
}

## forward simulations

H3N2.simulated.start.times <- matrix(data=0, ncol=length(H3N2_donor_names), nrow=n.sims)
colnames(H3N2.simulated.start.times) <- H3N2_donor_names

for (ferret in H3N2_donor_names){
  data <- H3N2_transmission_probs[[ferret]]
  bins <- round(seq(from=data[[1,1]], to=data[[length(data$dpe),1]], by=bin.length), digits=3)
  lambdas <- data %>%
    filter(dpe %in% bins) %>%
    ## no transmission if titer is at LOD
    filter(nw_titer > LOD) %>%
    pull(nw_titer) * H3N2.MLE
  transmission.probs <- 1 - exp(-lambdas*bin.length)
  ## stochastic simulations
  for (n in 1:n.sims){
    random.draws <- runif(length(transmission.probs), 0, 1)
    ## earliest simulated start time
    H3N2.simulated.start.times[n, ferret] <- bins[which(transmission.probs > random.draws)[1]]
  }
}

H3N2.simulated.start.times <- as.data.frame(H3N2.simulated.start.times)
H3N2.simulated.start.times$sim <- 1:n.sims
H3N2.simulated.start.times <- H3N2.simulated.start.times %>%
  pivot_longer(cols=1:length(H3N2_donor_names), names_to="animal", values_to="start.time")

## find first positive test for all contact ferrets
H3N2.first.positive.test <- data.frame()

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  H3N2.first.positive.test <- rbind(H3N2.first.positive.test, ferret_data[which.max(ferret_data$nw_titer > LOD),c("DI_RC_Pair","dpe")])
}
names(H3N2.first.positive.test) <- c("animal", "start.time")
H3N2.first.positive.test$animal <- as.character(H3N2.first.positive.test$animal)

panel_d <- ggplot() +
  geom_violin(data=H3N2.simulated.start.times, aes(x=animal, y=start.time)) +
  geom_point(data=H3N2.first.positive.test, aes(x=animal, y=start.time), color=H3N2_color) +
  labs(x="Index", y="Infection start time") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2)) +
  theme_light()

ggarrange(panel_b, panel_d, ncol=2, align="h")
