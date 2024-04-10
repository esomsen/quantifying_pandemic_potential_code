library(ggpubr)
library(khroma)


# H1N1 analysis -----------------------------------------------------------

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[1]

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

LOD <- 0.5

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

panel_a <- ggplot(H1N1_peak_titer, aes(x=donor_dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_peak_titer_regression)[[2]], 
              intercept = coef(H1N1_peak_titer_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  guides(color = "none") +
  xlab(NULL) +
  ylab(expression(paste("Peak titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## linear regresssion for time to peak titer
H1N1_time_peak$numeric_dose <- as.numeric(substr(H1N1_time_peak$donor_dose, 4, 4))
H1N1_time_peak_regression <- lm(time ~ numeric_dose, H1N1_time_peak)

panel_b <- ggplot(H1N1_time_peak, aes(x=donor_dose, y=time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_time_peak_regression)[[2]], 
              intercept = coef(H1N1_time_peak_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  guides(color = "none") +
  xlab(NULL) +
  ylab("Time to peak titer (days)") +
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
## linear regresssion for duration of infection
H1N1_duration$numeric_dose <- as.numeric(substr(H1N1_duration$donor_dose, 4, 4))
H1N1_duration_regression <- lm(duration ~ numeric_dose, H1N1_duration)

panel_c <- ggplot(H1N1_duration, aes(x=donor_dose, y=duration)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_duration_regression)[[2]], 
              intercept = coef(H1N1_duration_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  xlab(NULL) +
  ylab("Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## AUC - H1N1

## normalize titers by subtracting LOD and then calculate AUC

H1N1_AUCs <- c(rep(0, length(H1N1_recipient_names)))
names(H1N1_AUCs) <- H1N1_recipient_names
for (ferret in H1N1_recipient_names){
  ## filter data for individual ferrets and normalize by LOD
  ferret_data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    mutate(nw_titer = nw_titer - LOD)
  y_vals <- ferret_data$nw_titer
  x_vals <- ferret_data$dpe
  H1N1_AUCs[ferret] <- AUC(x=x_vals, y=y_vals, method = "trapezoid")
}

H1N1_AUC_df <- data.frame(Ferret_ID = names(H1N1_AUCs),
                          donor_dose = H1N1_time_peak$donor_dose,
                          AUC = H1N1_AUCs)

## linear regresssion for AUC
H1N1_AUC_df$numeric_dose <- as.numeric(substr(H1N1_AUC_df$donor_dose, 4, 4))
H1N1_AUC_regression <- lm(AUC ~ numeric_dose, H1N1_AUC_df)

panel_d <- ggplot(H1N1_AUC_df, aes(x=donor_dose, y=AUC)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_AUC_regression)[[2]], 
              intercept = coef(H1N1_AUC_regression)[[1]], 
              color=H1N1_color, linewidth=1) +
  ## add signif star
  annotate("text", label = "*", x=4, y=29, size=10) +
  xlab(NULL) +
  ylab("Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 4)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))


# H3N2 analysis -----------------------------------------------------------

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[2]

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

panel_e <- ggplot(H3N2_peak_titer, aes(x=donor_dose, y=nw_titer)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_peak_titer_regression)[[2]], 
              intercept = coef(H3N2_peak_titer_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  ## add signif star
  annotate("text", label = "*", x=4, y=7, size=10) +
  guides(color = "none") +
  xlab(expression(paste("Donor inoculum dose (", log[10], TCID[50], ")"))) +
  ylab(expression(paste("Peak titer (", log[10], TCID[50], ")"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8)) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## linear regresssion for time to peak titer
H3N2_time_peak$numeric_dose <- as.numeric(substr(H3N2_time_peak$donor_dose, 4, 4))
H3N2_time_peak_regression <- lm(time ~ numeric_dose, H3N2_time_peak)

panel_f <- ggplot(H3N2_time_peak, aes(x=donor_dose, y=time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_time_peak_regression)[[2]], 
              intercept = coef(H3N2_time_peak_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  guides(color = "none") +
  xlab(expression(paste("Donor inoculum dose (", log[10], TCID[50], ")"))) +
  ylab("Time to peak titer (days)") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10), limits = c(0, 10)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## duration of infection - H3N2

H3N2_infx_lengths <- c(rep(0, length(H3N2_recipient_names)))
names(H3N2_infx_lengths) <- H3N2_recipient_names

## treat all viral titers as "real"
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
  ## two ferrets have initial positive tests followed by negative tests
  ## deal with those separately
  if (ferret %in% c("2124", "626")){
    post_infection <- ferret_data %>%
      filter(dpe > 5)
    resolution_time <- post_infection[[which.max(post_infection$nw_titer == LOD), "dpe"]]
    duration <- resolution_time - first_positive_time 
  } else {
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
  }
  H3N2_infx_lengths[ferret] <- duration
}


H3N2_duration <- data.frame(Ferret_ID = H3N2_recipient_names, 
                            donor_dose = H3N2_time_peak$donor_dose, 
                            duration = H3N2_infx_lengths)
## linear regresssion for duration of infection
H3N2_duration$numeric_dose <- as.numeric(substr(H3N2_duration$donor_dose, 4, 4))
H3N2_duration_regression <- lm(duration ~ numeric_dose, H3N2_duration)

panel_g <- ggplot(H3N2_duration, aes(x=donor_dose, y=duration)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_duration_regression)[[2]], 
              intercept = coef(H3N2_duration_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  xlab(expression(paste("Donor inoculum dose (", log[10], TCID[50], ")"))) +
  ylab("Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6))

## AUC - H3N2

## normalize titers by subtracting LOD and then calculate AUC

H3N2_AUCs <- c(rep(0, length(H3N2_recipient_names)))
names(H3N2_AUCs) <- H3N2_recipient_names
for (ferret in H3N2_recipient_names){
  ## filter data for individual ferrets and normalize by LOD
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    mutate(nw_titer = nw_titer - LOD)
  y_vals <- ferret_data$nw_titer
  x_vals <- ferret_data$dpe
  H3N2_AUCs[ferret] <- AUC(x=x_vals, y=y_vals, method = "trapezoid")
}

H3N2_AUC_df <- data.frame(Ferret_ID = names(H3N2_AUCs),
                          donor_dose = H3N2_time_peak$donor_dose,
                          AUC = H3N2_AUCs)

## linear regresssion for AUC
H3N2_AUC_df$numeric_dose <- as.numeric(substr(H3N2_AUC_df$donor_dose, 4, 4))
H3N2_AUC_regression <- lm(AUC ~ numeric_dose, H3N2_AUC_df)

panel_h <- ggplot(H3N2_AUC_df, aes(x=donor_dose, y=AUC)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_AUC_regression)[[2]], 
              intercept = coef(H3N2_AUC_regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  xlab(expression(paste("Donor inoculum dose (", log[10], TCID[50], ")"))) +
  ylab("Duration of infection (days)") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 4)) +
  theme_light() +
  scale_x_discrete(limits=c("10^0","10^1","10^2", "10^3", "10^4", "10^5", "10^6"),
                   labels=c(0, 1, 2, 3, 4, 5, 6)) 

# full figure -------------------------------------------------------------

col1 <- ggarrange(panel_a, panel_e, ncol=1, nrow=2, labels = c("A", "E"))
col2 <- ggarrange(panel_b, panel_f, ncol=1, nrow=2, labels = c("B", "F"))
col3 <- ggarrange(panel_c, panel_g, ncol=1, nrow=2, labels = c("C", "G"))
col4 <- ggarrange(panel_d, panel_h, ncol=1, nrow=2, labels = c("D", "H"))

ggarrange(col1, col2, col3, col4, nrow=1, ncol=4)
