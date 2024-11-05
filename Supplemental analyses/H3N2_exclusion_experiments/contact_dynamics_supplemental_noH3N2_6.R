library(khroma)
library(tidyverse)
library(ggpubr)
library(DescTools)

# H3N2 analyses -----------------------------------------------------------

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[[2]]

LOD <- 0.5

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "dose", "dpi", "nw_titer", "donor_dose")

H3N2_RC_ferrets <- H3N2_ferrets %>%
  filter(DI_RC == "RC") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, donor_dose, dpi, nw_titer) %>%
  ## exclude all negative tests at 13 and 14 dpi
  filter(dpi < 14) %>%
  mutate(dpe = dpi - 1) %>%
  ## remove 10^6 ferrets
  filter(donor_dose != "10^6")

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
  mutate(dpe = dpi - 1) %>%
  ## remove 10^6 ferrets
  filter(dose != "10^6")

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

## stays insignificant

panel_a <- ggplot(H3N2.AUC.AUC, aes(x=Index.AUC, y=Contact.AUC)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  geom_abline(slope = coef(H3N2_DI_RC_AUC_regression)[[2]], 
              intercept = coef(H3N2_DI_RC_AUC_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  scale_y_continuous(breaks=seq(0, 30, 5), limits = c(0, 30)) +
  scale_x_continuous(breaks=seq(0, 30, 5), limits = c(0, 30)) +
  labs(title="A", x="Index AUC", y="Contact AUC") +
  theme_light()

## calculate donor AUC

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

## peak titer timing and magnitude - H3N2

H3N2_peak_titer <- data.frame()
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
  H3N2_peak_titer <- rbind(H3N2_peak_titer, ferret_data[max_titer_index,])
}

H3N2_time_peak <- merge(H3N2_time_peak, H3N2_donor_AUCs, by="Ferret_ID")
H3N2_peak_titer <- merge(H3N2_peak_titer[,c(1,2,4)], H3N2_donor_AUCs, by="Ferret_ID")

## linear regression for peak titer value
## use donor AUC
H3N2_peak_titer_regression <- lm(nw_titer ~ donor_AUC, H3N2_peak_titer)

## becomes insignificant at a=0.05

panel_b <- ggplot(H3N2_peak_titer, aes(x=donor_AUC, y=nw_titer)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_peak_titer_regression)[[2]], 
              intercept = coef(H3N2_peak_titer_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  guides(color = "none") +
  labs(title="B", x="Index AUC", y=expression(paste("Peak titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  scale_x_continuous(breaks=seq(0, 30, 5), limits = c(0, 30)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2)

## linear regression for time to peak titer
H3N2_time_peak_regression <- lm(time ~ donor_AUC, H3N2_time_peak)

## stays insignificant 

panel_c <- ggplot(H3N2_time_peak, aes(x=donor_AUC, y=time)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_time_peak_regression)[[2]], 
              intercept = coef(H3N2_time_peak_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  guides(color = "none") +
  labs(title="C", x="Index AUC", y="Time to peak titer (days)") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10), limits = c(0, 10)) +
  scale_x_continuous(breaks=seq(0, 30, 5), limits = c(0, 30)) +
  theme_light()

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

## stays insignificant

panel_d <- ggplot(H3N2_infx_lengths, aes(x=donor_AUC, y=duration)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_duration_regression)[[2]], 
              intercept = coef(H3N2_duration_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  labs(title = "D", x="Index AUC", y="Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks=seq(0, 30, 5), limits = c(0, 30)) +
  theme_light()

ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, align="h")
