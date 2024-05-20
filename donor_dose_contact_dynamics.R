library(ggpubr)
library(khroma)
library(tidyverse)
library(DescTools)

# H1N1 analysis -----------------------------------------------------------

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[[1]]

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

## initial titer (0dpe) of index animals and donor AUC

H1N1.init.titers <- data.frame()
H1N1.donor.AUC <- data.frame()

for (ferret in H1N1_donor_names){
  ferret_data <- H1N1_DI_ferrets %>%
    filter(Ferret_ID == ferret)
  H1N1.init.titers <- rbind(H1N1.init.titers, ferret_data[1,])
  ## normalize titers by subtracting LOD and then calculate AUC
  y_vals <- ferret_data$nw_titer - LOD
  x_vals <- ferret_data$dpe
  H1N1.donor.AUC <- rbind(H1N1.donor.AUC, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## linear regression for initial titers
H1N1.init.titers$numeric_dose <- as.numeric(substr(H1N1.init.titers$dose, 4, 4))
H1N1.init.regression <- lm(nw_titer ~ numeric_dose, H1N1.init.titers)

panel_f <- ggplot(H1N1.init.titers, aes(x=dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill="black", color="black") +
  ## add regression line
  geom_abline(slope = coef(H1N1.init.regression)[[2]], 
              intercept = coef(H1N1.init.regression)[[1]], 
              color="black", linewidth=1) +
  guides(color = "none") +
  labs(title="F", xlab=NULL, y=expression(paste("Initial titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## linear regression for index AUC

H1N1.donor.AUC$numeric_dose <- as.numeric(substr(H1N1.donor.AUC$dose, 4, 4))
H1N1.donor.AUC.regression <- lm(AUC ~ numeric_dose, H1N1.donor.AUC)

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

## plot index and contact together
all.AUC <- data.frame("dose" = c(H1N1.donor.AUC[,"dose"], H1N1.AUCs[,"donor_dose"]), 
                      "AUC" = c(H1N1.donor.AUC[,"AUC"], H1N1.AUCs[,"AUC"]))
#all.AUC <- rbind(H1N1.donor.AUC[,c("dose", "AUC")], H1N1.AUCs[,c("_dose", "AUC")])
all.AUC$type <- c(rep("donor", length(H1N1_donor_names)), rep("recipient", length(H1N1_recipient_names)))

panel_g <- ggplot(all.AUC, aes(x=dose, y=AUC, color=type, fill=type)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  scale_color_manual(values = c(donor="black", recipient=H1N1_color)) +
  scale_fill_manual(name=NULL, labels = c("Index", "Contact"), values = c("black", H1N1_color)) +
  ## add regression line
  geom_abline(slope = coef(H1N1_AUC_regression)[[2]], 
              intercept = coef(H1N1_AUC_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  geom_abline(slope = coef(H1N1.donor.AUC.regression)[[2]], 
              intercept = coef(H1N1.donor.AUC.regression)[[1]], 
              color="black", linewidth=1) +
  labs(title="G", x=NULL, y="AUC (log-transformed)") +
  guides(color="none") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 4)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## time of initial positive test in contacts

H1N1.time.positive <- data.frame()

for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## find first positive test
  first_positive_time <- ferret_data[[which.max(ferret_data$nw_titer > LOD), "dpe"]]
  H1N1.time.positive <- rbind(H1N1.time.positive, c(ferret_data[1,1:2], "pos_time" = first_positive_time))
}

## linear regression for time to first positive test
H1N1.time.positive$numeric_dose <- as.numeric(substr(H1N1.time.positive$donor_dose, 4, 4))
H1N1.time.positive.regression <- lm(pos_time ~ numeric_dose, H1N1.time.positive)

panel_h <- ggplot(H1N1.time.positive, aes(x=donor_dose, y=pos_time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1.time.positive.regression)[[2]], 
              intercept = coef(H1N1.time.positive.regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  labs(title="H", x=NULL, y="Time to first positive test (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## peak titer timing and magnitude - H1N1

H1N1_peak_titer <- data.frame()
H1N1_time_peak <- data.frame()

for (ferret in H1N1_recipient_names){
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  max_titer_index <- which.max(ferret_data$nw_titer)
  first_positive_index <- which.max(ferret_data$nw_titer > LOD)
  ## time since last negative test to the peak test titer
  if (first_positive_index > 1){
    time <- ferret_data[[max_titer_index,"dpe"]] - ferret_data[[first_positive_index-1,"dpe"]]
  } else { ## if the first positive test is the first test, consider t=0 to be start of infection
    time <- ferret_data[[max_titer_index,"dpe"]]
  }
  H1N1_time_peak <- rbind(H1N1_time_peak, cbind(ferret_data[1, 1:2], time))
  H1N1_peak_titer <- rbind(H1N1_peak_titer, ferret_data[max_titer_index,])
}

## linear regression for peak titer value
## use donor inoculum dose as a numeric value on log-scale
H1N1_peak_titer$numeric_dose <- as.numeric(substr(H1N1_peak_titer$donor_dose, 4, 4))
H1N1_peak_titer_regression <- lm(nw_titer ~ numeric_dose, H1N1_peak_titer)

panel_i <- ggplot(H1N1_peak_titer, aes(x=donor_dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_peak_titer_regression)[[2]], 
              intercept = coef(H1N1_peak_titer_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  guides(color = "none") +
  labs(title="I", x=NULL, y=expression(paste("Peak titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## linear regression for time to peak titer
H1N1_time_peak$numeric_dose <- as.numeric(substr(H1N1_time_peak$donor_dose, 4, 4))
H1N1_time_peak_regression <- lm(time ~ numeric_dose, H1N1_time_peak)

panel_j <- ggplot(H1N1_time_peak, aes(x=donor_dose, y=time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_time_peak_regression)[[2]], 
              intercept = coef(H1N1_time_peak_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  guides(color = "none") +
  labs(title="J", x=NULL, y="Time to peak titer (days)") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10), limits = c(0, 10)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

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
                            donor_dose = H1N1_time_peak$donor_dose, 
                            duration = H1N1_infx_lengths)
## linear regression for duration of infection
H1N1_duration$numeric_dose <- as.numeric(substr(H1N1_duration$donor_dose, 4, 4))
H1N1_duration_regression <- lm(duration ~ numeric_dose, H1N1_duration)

panel_k <- ggplot(H1N1_duration, aes(x=donor_dose, y=duration)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_duration_regression)[[2]], 
              intercept = coef(H1N1_duration_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  labs(title = "K", x=NULL, y="Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

H1N1_kinetcs_plot <- ggarrange(panel_f, panel_g, panel_h, panel_i, panel_j, panel_k, nrow=2, ncol=3, align = "v", common.legend = T, legend="right")
H1N1_kinetcs_plot <- annotate_figure(H1N1_kinetcs_plot, bottom = text_grob(expression(paste("Inoculum dose (", log[10], TCID[50], ")"))))

# H3N2 analysis -----------------------------------------------------------

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

## initial titer (0dpe) of index animals and donor AUC

H3N2.init.titers <- data.frame()
H3N2.donor.AUC <- data.frame()

for (ferret in H3N2_donor_names){
  ferret_data <- H3N2_DI_ferrets %>%
    filter(Ferret_ID == ferret)
  H3N2.init.titers <- rbind(H3N2.init.titers, ferret_data[1,])
  ## normalize titers by subtracting LOD and then calculate AUC
  y_vals <- ferret_data$nw_titer - LOD
  x_vals <- ferret_data$dpe
  H3N2.donor.AUC <- rbind(H3N2.donor.AUC, c(ferret_data[1,], "AUC"=AUC(x=x_vals, y=y_vals, method = "trapezoid")))
}

## linear regression for initial titers
H3N2.init.titers$numeric_dose <- as.numeric(substr(H3N2.init.titers$dose, 4, 4))
H3N2.init.regression <- lm(nw_titer ~ numeric_dose, H3N2.init.titers)

panel_g <- ggplot(H3N2.init.titers, aes(x=dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill="black", color="black") +
  ## add regression line
  geom_abline(slope = coef(H3N2.init.regression)[[2]], 
              intercept = coef(H3N2.init.regression)[[1]], 
              color="black", linewidth=1) +
  guides(color = "none") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6)) +
  labs(title="G", x=NULL, y=expression(paste("Initial titer (", log[10], TCID[50], ")")))

## linear regression for index AUC

H3N2.donor.AUC$numeric_dose <- as.numeric(substr(H3N2.donor.AUC$dose, 4, 4))
H3N2.donor.AUC.regression <- lm(AUC ~ numeric_dose, H3N2.donor.AUC)

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

## plot index and contact together
all.AUC <- data.frame("dose" = c(H3N2.donor.AUC[,"dose"], H3N2.AUCs[,"donor_dose"]), 
                      "AUC" = c(H3N2.donor.AUC[,"AUC"], H3N2.AUCs[,"AUC"]))
#all.AUC <- rbind(H3N2.donor.AUC[,c("dose", "AUC")], H3N2.AUCs[,c("_dose", "AUC")])
all.AUC$type <- c(rep("donor", length(H3N2_donor_names)), rep("recipient", length(H3N2_recipient_names)))

panel_h <- ggplot(all.AUC, aes(x=dose, y=AUC, color=type, fill=type)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  scale_color_manual(values = c(donor="black", recipient=H3N2_color)) +
  scale_fill_manual(name=NULL, labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  ## add regression line
  geom_abline(slope = coef(H3N2_AUC_regression)[[2]], 
              intercept = coef(H3N2_AUC_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  geom_abline(slope = coef(H3N2.donor.AUC.regression)[[2]], 
              intercept = coef(H3N2.donor.AUC.regression)[[1]], 
              color="black", linewidth=1) +
  labs(title="H", x=NULL, y="AUC (log-transformed)") +
  guides(color="none") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 4)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

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

panel_i <- ggplot(H3N2.time.positive, aes(x=donor_dose, y=pos_time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2.time.positive.regression)[[2]], 
              intercept = coef(H3N2.time.positive.regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  labs(title="I", x=NULL, y="Time to first positive test (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

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

## linear regression for peak titer value
## use donor inoculum dose as a numeric value on log-scale
H3N2_peak_titer$numeric_dose <- as.numeric(substr(H3N2_peak_titer$donor_dose, 4, 4))
H3N2_peak_titer_regression <- lm(nw_titer ~ numeric_dose, H3N2_peak_titer)

panel_j <- ggplot(H3N2_peak_titer, aes(x=donor_dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_peak_titer_regression)[[2]], 
              intercept = coef(H3N2_peak_titer_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  guides(color = "none") +
  labs(title="J", x=NULL, y=expression(paste("Peak titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## linear regression for time to peak titer
H3N2_time_peak$numeric_dose <- as.numeric(substr(H3N2_time_peak$donor_dose, 4, 4))
H3N2_time_peak_regression <- lm(time ~ numeric_dose, H3N2_time_peak)

panel_k <- ggplot(H3N2_time_peak, aes(x=donor_dose, y=time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_time_peak_regression)[[2]], 
              intercept = coef(H3N2_time_peak_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  guides(color = "none") +
  labs(title="K", x=NULL, y="Time to peak titer (days)") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10), limits = c(0, 10)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

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
    ## provided that there isn't an initial positive test (happens at 3, 5dpe)
    post_infection <- post_infection %>%
      filter(dpe > 5)
    resolution_time <- post_infection[[which.max(post_infection$nw_titer == LOD), "dpe"]]
  }
  duration <- resolution_time - first_positive_time
  H3N2_infx_lengths[ferret] <- duration
}

H3N2_duration <- data.frame(Ferret_ID = H3N2_recipient_names, 
                            donor_dose = H3N2_time_peak$donor_dose, 
                            duration = H3N2_infx_lengths)
## linear regression for duration of infection
H3N2_duration$numeric_dose <- as.numeric(substr(H3N2_duration$donor_dose, 4, 4))
H3N2_duration_regression <- lm(duration ~ numeric_dose, H3N2_duration)

panel_l <- ggplot(H3N2_duration, aes(x=donor_dose, y=duration)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_duration_regression)[[2]], 
              intercept = coef(H3N2_duration_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  labs(title = "L", x=NULL, y="Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

H3N2_kinetcs_plot <- ggarrange(panel_g, panel_h, panel_i, panel_j, panel_k, panel_l, nrow=2, ncol=3, align = "v", common.legend = T, legend="right")
H3N2_kinetcs_plot <- annotate_figure(H3N2_kinetcs_plot, bottom = text_grob(expression(paste("Inoculum dose (", log[10], TCID[50], ")"))))

