library(khroma)
library(tidyverse)
library(DescTools)

H1N1_ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H1N1_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

H1N1_RC_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer) %>%
  ## removing the uninfected ferret
  filter(Ferret_ID != 7824) %>%
  mutate(dpe = dpi - 1)
H1N1_recipient_names <- unique(H1N1_RC_ferrets$Ferret_ID)

LOD <- 0.5

# H1N1 --------------------------------------------------------------------

## ferret 5775; remove "rebound" positive test at 11 dpe

H1N1_RC_ferrets[69, 4] <- LOD

## AUC - H1N1 contacts

H1N1.AUCs <- data.frame()

for (ferret in H1N1_recipient_names){
  ## filter data for individual ferrets and normalize by LOD
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    ## normalize titers by subtracting LOD and then calculate AUC
    mutate(nw_titer = nw_titer - LOD)
  y_vals <- ferret_data$nw_titer
  x_vals <- ferret_data$dpe
  H1N1.AUCs <- rbind(H1N1.AUCs, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## linear regresssion for recipient AUC
H1N1.AUCs$numeric_dose <- as.numeric(substr(H1N1.AUCs$donor_dose, 4, 4))
H1N1_AUC_regression <- lm(AUC ~ numeric_dose, H1N1.AUCs)

## AUC stays significant, slope changes by 0.02

## duration of infection - H1N1

H1N1_infx_lengths <- c(rep(0, length(H1N1_recipient_names)))
names(H1N1_infx_lengths) <- H1N1_recipient_names

for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## find first positive test
  first_positive_index <- which.max(ferret_data$nw_titer > LOD)
  ## assume that viral growth begins immediately following the last negative test
  if (first_positive_index > 1){
    first_positive_time <- ferret_data[[first_positive_index-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    first_positive_time <- 0
  }
  ## find the time infection resolves
  post_infection <- ferret_data %>%
    filter(dpe > first_positive_time)
  if (post_infection[[length(post_infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    resolution_time <- post_infection[[length(post_infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    resolution_time <- post_infection[[which.max(post_infection$nw_titer == LOD), "dpe"]]
  }
  duration <- resolution_time - first_positive_time
  H1N1_infx_lengths[ferret] <- duration
}

H1N1_duration <- data.frame(Ferret_ID = H1N1_recipient_names, 
                            donor_dose = H1N1.AUCs$donor_dose, 
                            duration = H1N1_infx_lengths)
## linear regression for duration of infection
H1N1_duration$numeric_dose <- as.numeric(substr(H1N1_duration$donor_dose, 4, 4))
H1N1_duration_regression <- lm(duration ~ numeric_dose, H1N1_duration)

## duration of infection remains insignificant, slope changes by 0.03

# H3N2 --------------------------------------------------------------------

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

LOD <- 0.5

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


## ferret 626; remove initial positive test followed by a negative test

H3N2_RC_ferrets[49, 4] <- LOD

## ferret 2124; remove initial positive test followed by a negative test

H3N2_RC_ferrets[2,4] <- LOD

## AUC - H3N2 contacts

H3N2.AUCs <- data.frame()

for (ferret in H3N2_recipient_names){
  ## filter data for individual ferrets and normalize by LOD
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    ## normalize titers by subtracting LOD and then calculate AUC
    mutate(nw_titer = nw_titer - LOD)
  y_vals <- ferret_data$nw_titer
  x_vals <- ferret_data$dpe
  H3N2.AUCs <- rbind(H3N2.AUCs, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## linear regresssion for recipient AUC
H3N2.AUCs$numeric_dose <- as.numeric(substr(H3N2.AUCs$donor_dose, 4, 4))
H3N2_AUC_regression <- lm(AUC ~ numeric_dose, H3N2.AUCs)

## AUC closer to significant, p=0.0857

## time of initial positive test in contacts

H3N2.time.positive <- data.frame()

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## find first positive test
  first_positive_time <- ferret_data[[which.max(ferret_data$nw_titer > LOD), "dpe"]]
  H3N2.time.positive <- rbind(H3N2.time.positive, c(ferret_data[1,1:2], "pos_time" = first_positive_time))
}

## linear regression for time to first positive test
H3N2.time.positive$numeric_dose <- as.numeric(substr(H3N2.time.positive$donor_dose, 4, 4))
H3N2.time.positive.regression <- lm(pos_time ~ numeric_dose, H3N2.time.positive)

## time to first positive test remains insignificant, slope changes by 0.01

## duration of infection - H3N2

H3N2_infx_lengths <- c(rep(0, length(H3N2_recipient_names)))
names(H3N2_infx_lengths) <- H3N2_recipient_names

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## find first positive test
  first_positive_index <- which.max(ferret_data$nw_titer > LOD)
  ## assume that viral growth begins immediately following the last negative test
  if (first_positive_index > 1){
    first_positive_time <- ferret_data[[first_positive_index-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    first_positive_time <- 0
  }
  ## find the time infection resolves
  post_infection <- ferret_data %>%
    filter(dpe > first_positive_time)
  if (post_infection[[length(post_infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    resolution_time <- post_infection[[length(post_infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    resolution_time <- post_infection[[which.max(post_infection$nw_titer == LOD), "dpe"]]
  }
  duration <- resolution_time - first_positive_time
  H3N2_infx_lengths[ferret] <- duration
}

H3N2_duration <- data.frame(Ferret_ID = H3N2_recipient_names, 
                            donor_dose = H3N2.time.positive$donor_dose, 
                            duration = H3N2_infx_lengths)
## linear regression for duration of infection
H3N2_duration$numeric_dose <- as.numeric(substr(H3N2_duration$donor_dose, 4, 4))
H3N2_duration_regression <- lm(duration ~ numeric_dose, H3N2_duration)

## duration remains insignificant
