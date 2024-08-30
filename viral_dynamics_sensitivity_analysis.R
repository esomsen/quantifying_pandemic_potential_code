library(khroma)
library(tidyverse)
library(DescTools)
library(ggpubr)

# H1N1 --------------------------------------------------------------------

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

H1N1_DI_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "DI") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, dose, dpi, nw_titer) %>%
  mutate(dpe = dpi - 1)
H1N1_donor_names <- unique(H1N1_DI_ferrets$Ferret_ID)

LOD <- 0.5

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[[1]]

## ferret 5775; remove "rebound" positive test at 11 dpe

## titers figure

donor_5774 <- H1N1_DI_ferrets %>%
  filter(Ferret_ID == "5774")
recipient_5775 <- H1N1_RC_ferrets %>%
  filter(Ferret_ID == "5775")

panel_a_titers <- data.frame(donor = donor_5774$nw_titer, 
                             recipient = recipient_5775$nw_titer) %>%
  pivot_longer(cols=1:2, names_to="animal", values_to="titer") %>%
  mutate(dpe = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

panel_a <- ggplot(panel_a_titers, aes(x=dpe, y=titer, color=animal)) +
  geom_point(shape=15, size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H1N1_color)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  geom_segment(aes(x=7, y=4, xend=8.5, yend=2), arrow = arrow(length=unit(.5, 'cm')), color="black", lwd=2) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  labs(x = "Days post exposure", y=expression(paste("Viral titer (", log[10], TCID[50], ")")), color=NULL)

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

## AUC - H1N1 donors

H1N1.donor.AUC <- data.frame()

for (ferret in H1N1_donor_names){
  ferret_data <- H1N1_DI_ferrets %>%
    filter(Ferret_ID == ferret)
  ## normalize titers by subtracting LOD and then calculate AUC
  y_vals <- ferret_data$nw_titer - LOD
  x_vals <- ferret_data$dpe
  H1N1.donor.AUC <- rbind(H1N1.donor.AUC, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## AUC vs AUC

H1N1.AUC.AUC <- merge(x=H1N1.donor.AUC[,c(1, 6)], y=H1N1_ferrets[,c(1,3,4)], by="Ferret_ID", sort=F, all.x=T) %>%
  unique()
H1N1.AUC.AUC <- merge(x=H1N1.AUC.AUC, y=H1N1.AUCs[,c(1,6)], by.x="DI_RC_Pair", by.y="Ferret_ID", sort=F)
names(H1N1.AUC.AUC) <- c("Contact", "Index", "Index.AUC", "donor_dose", "Contact.AUC")
H1N1_DI_RC_AUC_regression <- lm(Contact.AUC ~ Index.AUC, H1N1.AUC.AUC)

## AUC remains insignificant

## calculate donor AUC

H1N1_donor_AUCs <- data.frame()

for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  first_positive_time <- ferret_data[[which.max(ferret_data$nw_titer > LOD), "dpe"]]
  ## calculate AUC of donor up to first positive test
  donor_name <- H1N1_ferrets[[which.max(H1N1_ferrets$DI_RC_Pair == ferret),"Ferret_ID"]]
  donor_data <- H1N1_DI_ferrets %>%
    filter(Ferret_ID == donor_name)
  y_vals <- donor_data$nw_titer - LOD
  x_vals <- donor_data$dpe
  H1N1_donor_AUCs <- rbind(H1N1_donor_AUCs, c(ferret_data[1,1], "donor_AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## duration of infection - H1N1

H1N1_infx_lengths <- data.frame()

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
  H1N1_infx_lengths <- rbind(H1N1_infx_lengths, cbind(ferret_data[1, 1:2], "duration"=resolution_time - first_positive_time))
}

H1N1_infx_lengths <- merge(H1N1_infx_lengths, H1N1_donor_AUCs, by="Ferret_ID")

## linear regression for duration of infection
H1N1_duration_regression <- lm(duration ~ donor_AUC, H1N1_infx_lengths)

## duration of infection remains insignificant, slope changes by 0.03

# H3N2 --------------------------------------------------------------------

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

LOD <- 0.5

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[[2]]

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

## ferret 626; remove initial positive test followed by a negative test

donor_671 <- H3N2_DI_ferrets %>%
  filter(Ferret_ID == "671")
recipient_626 <- H3N2_RC_ferrets %>%
  filter(Ferret_ID == "626")

panel_b_titers <- data.frame(donor = donor_671$nw_titer[1:6], 
                             recipient = recipient_626$nw_titer) %>%
  pivot_longer(cols=1:2, names_to="animal", values_to="titer") %>%
  mutate(dpe = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) %>%
  mutate(LOD_shape = ifelse(titer <= 0.5, "below", "above"))

panel_b <- ggplot(panel_b_titers, aes(x=dpe, y=titer, color=animal)) +
  geom_point(shape=16, size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  geom_segment(aes(x=3, y=4, xend=2, yend=2), arrow = arrow(length=unit(.5, 'cm')), color="black", lwd=2) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  labs(x = "Days post exposure", y=expression(paste("Viral titer (", log[10], TCID[50], ")")), color=NULL)

H3N2_RC_ferrets[49, 4] <- LOD

## ferret 2124; remove initial positive test followed by a negative test

donor_2120 <- H3N2_DI_ferrets %>%
  filter(Ferret_ID == "2120")
recipient_2124 <- H3N2_RC_ferrets %>%
  filter(Ferret_ID == "2124")

panel_c_titers <- data.frame(donor = donor_2120$nw_titer, 
                             recipient = recipient_2124$nw_titer) %>%
  pivot_longer(cols=1:2, names_to="animal", values_to="titer") %>%
  mutate(dpe = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) %>%
  mutate(LOD_shape = ifelse(titer <= 0.5, "below", "above"))

panel_c <- ggplot(panel_c_titers, aes(x=dpe, y=titer, color=animal)) +
  geom_point(shape=16, size=2) +
  geom_line(linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  geom_segment(aes(x=4, y=3, xend=2, yend=1.5), arrow = arrow(length=unit(.5, 'cm')), color="black", lwd=2) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  labs(x = "Days post exposure", y=expression(paste("Viral titer (", log[10], TCID[50], ")")), color=NULL)


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

## AUC - H3N2 donors

H3N2.donor.AUC <- data.frame()

for (ferret in H3N2_donor_names){
  ferret_data <- H3N2_DI_ferrets %>%
    filter(Ferret_ID == ferret)
  ## normalize titers by subtracting LOD and then calculate AUC
  y_vals <- ferret_data$nw_titer - LOD
  x_vals <- ferret_data$dpe
  H3N2.donor.AUC <- rbind(H3N2.donor.AUC, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## AUC vs AUC

H3N2.AUC.AUC <- merge(x=H3N2.donor.AUC[,c(1, 6)], y=H3N2_ferrets[,c(1,3,4)], by="Ferret_ID", sort=F, all.x=T) %>%
  unique()
H3N2.AUC.AUC <- merge(x=H3N2.AUC.AUC, y=H3N2.AUCs[,c(1,6)], by.x="DI_RC_Pair", by.y="Ferret_ID", sort=F)
names(H3N2.AUC.AUC) <- c("Contact", "Index", "Index.AUC", "donor_dose", "Contact.AUC")
H3N2_DI_RC_AUC_regression <- lm(Contact.AUC ~ Index.AUC, H3N2.AUC.AUC)

## AUC closer to significant, p=0.0857

## calculate donor AUC (up to time of first positive test?)

H3N2_donor_AUCs <- data.frame()

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  first_positive_time <- ferret_data[[which.max(ferret_data$nw_titer > LOD), "dpe"]]
  ## calculate AUC of donor up to first positive test
  donor_name <- H3N2_ferrets[[which.max(H3N2_ferrets$DI_RC_Pair == ferret),"Ferret_ID"]]
  donor_data <- H3N2_DI_ferrets %>%
    filter(Ferret_ID == donor_name)
  y_vals <- donor_data$nw_titer - LOD
  x_vals <- donor_data$dpe
  H3N2_donor_AUCs <- rbind(H3N2_donor_AUCs, c(ferret_data[1,1], "donor_AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## time of initial positive test in contacts

H3N2_time_peak <- data.frame()

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  max_titer_index <- which.max(ferret_data$nw_titer)
  first_positive_index <- which.max(ferret_data$nw_titer > LOD)
  ## time since last negative test to the peak test titer
  if (first_positive_index > 1){
    time <- ferret_data[[max_titer_index,"dpe"]] - ferret_data[[first_positive_index-1,"dpe"]]
  } else { ## if the first positive test is the first test, consider t=0 to be start of infection
    time <- ferret_data[[max_titer_index,"dpe"]]
  }
  H3N2_time_peak <- rbind(H3N2_time_peak, cbind(ferret_data[1, 1:2], time))
}

H3N2_time_peak <- merge(H3N2_time_peak, H3N2_donor_AUCs, by="Ferret_ID")

## linear regression for time to peak titer
H3N2_time_peak_regression <- lm(time ~ donor_AUC, H3N2_time_peak)

## time to first positive test remains insignificant, slope changes by 0.01

## duration of infection - H3N2

H3N2_infx_lengths <- data.frame()

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
  H3N2_infx_lengths <- rbind(H3N2_infx_lengths, cbind(ferret_data[1, 1:2], "duration"=resolution_time - first_positive_time))
}

H3N2_infx_lengths <- merge(H3N2_infx_lengths, H3N2_donor_AUCs, by="Ferret_ID")

## linear regression for duration of infection
H3N2_duration_regression <- lm(duration ~ donor_AUC, H3N2_infx_lengths)

## duration remains insignificant


## combo plots

ggarrange(panel_a, panel_b, panel_c, legend="bottom", ncol=3, labels=c("A", "B", "C"))
