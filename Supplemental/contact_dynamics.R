library(khroma)
library(tidyverse)
library(ggpubr)
library(DescTools)

LOD <- 1

# H1N1 analyses -----------------------------------------------------------

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[[1]]

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

H1N1_DI_ferrets <- H1N1_ferrets %>%
  filter(DI_RC == "DI") %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  dplyr::select(Ferret_ID, Dose, dpi, nw_titer) %>%
  mutate(dpe = dpi - 1)
H1N1_donor_names <- unique(H1N1_DI_ferrets$Ferret_ID)

## peak titer H1N1

H1N1.peak.titer <- data.frame()

for (ferret in H1N1_recipient_names){
  tmp.data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  H1N1.peak.titer <- rbind(H1N1.peak.titer, tmp.data[which.max(tmp.data$nw_titer),])
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## linear regression for peak titer value
## make dose a numeric
H1N1.peak.titer$numeric_dose <- c(6, 6, 6, 6, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0)
H1N1_peak_titer_regression <- lm(nw_titer ~ numeric_dose, H1N1.peak.titer)

panel_a <- ggplot(H1N1.peak.titer, aes(x=donor_dose, y=nw_titer)) +
  geom_point(size=2, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_peak_titer_regression)[[2]], 
              intercept = coef(H1N1_peak_titer_regression)[[1]], 
              color="black", linewidth=1) +
  guides(color = "none") +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y=expression(paste("Peak viral titer (", TCID[50],"/mL", ")"))) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2)

## duration of infection - H1N1

H1N1.infx.lengths <- data.frame()

for (ferret in H1N1_recipient_names){
  tmp.data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## assume that viral growth begins immediately following the last negative test
  if (which.max(tmp.data$nw_titer > LOD) > 1){
    tmp.pos.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    tmp.pos.time <- 0
  }
  ## find the time infection resolves
  tmp.post.infection <- tmp.data %>%
    filter(dpe > tmp.pos.time)
  if (tmp.post.infection[[length(tmp.post.infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    tmp.end.time <- tmp.post.infection[[length(tmp.post.infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    tmp.end.time <- tmp.post.infection[[which.max(tmp.post.infection$nw_titer <= LOD), "dpe"]]
  }
  H1N1.infx.lengths <- rbind(H1N1.infx.lengths, cbind(tmp.data[1, 1:2], "duration"=tmp.end.time - tmp.pos.time))
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## linear regression for duration of infection
H1N1.infx.lengths$numeric_dose <- c(6, 6, 6, 6, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0)
H1N1_duration_regression <- lm(duration ~ numeric_dose, H1N1.infx.lengths)

panel_b <- ggplot(H1N1.infx.lengths, aes(x=donor_dose, y=duration)) +
  geom_point(size=2, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_duration_regression)[[2]], 
              intercept = coef(H1N1_duration_regression)[[1]], 
              color="black", linewidth=1) +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y="Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  theme_light()

# H3N2 analyses -----------------------------------------------------------

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[[2]]

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

## peak titer H3N2

H3N2.peak.titer <- data.frame()

for (ferret in H3N2_recipient_names){
  tmp.data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  H3N2.peak.titer <- rbind(H3N2.peak.titer, tmp.data[which.max(tmp.data$nw_titer),])
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## linear regression for peak titer value
H3N2.peak.titer$numeric_dose <- c(6, 6, 6, 4, 4, 4, 4, 3, 2, 2, 1)
H3N2_peak_titer_regression <- lm(nw_titer ~ numeric_dose, H3N2.peak.titer)

panel_c <- ggplot(H3N2.peak.titer, aes(x=donor_dose, y=nw_titer)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_peak_titer_regression)[[2]], 
             intercept = coef(H3N2_peak_titer_regression)[[1]], 
              color="black", linewidth=1) +
  guides(color = "none") +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y=expression(paste("Peak viral titer (", TCID[50],"/mL", ")"))) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2)

## duration of infection - H3N2

H3N2.infx.lengths <- data.frame()

for (ferret in H3N2_recipient_names){
  tmp.data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## assume that viral growth begins immediately following the last negative test
  if (which.max(tmp.data$nw_titer > LOD) > 1){
    tmp.pos.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    tmp.pos.time <- 0
  }
  ## find the time infection resolves
  tmp.post.infection <- tmp.data %>%
    filter(dpe > tmp.pos.time)
  if (tmp.post.infection[[length(tmp.post.infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    tmp.end.time <- tmp.post.infection[[length(tmp.post.infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    tmp.end.time <- tmp.post.infection[[which.max(tmp.post.infection$nw_titer <= LOD), "dpe"]]
  }
  H3N2.infx.lengths <- rbind(H3N2.infx.lengths, cbind(tmp.data[1, 1:2], "duration"=tmp.end.time - tmp.pos.time))
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}
## fix value for 2124
H3N2.infx.lengths[2,3] <- 11

## linear regression for duration of infection
H3N2.infx.lengths$numeric_dose <- c(6, 6, 6, 4, 4, 4, 4, 3, 2, 2, 1)
H3N2_duration_regression <- lm(duration ~ numeric_dose, H3N2.infx.lengths)

panel_d <- ggplot(H3N2.infx.lengths, aes(x=donor_dose, y=duration)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_duration_regression)[[2]], 
             intercept = coef(H3N2_duration_regression)[[1]], 
             color="black", linewidth=1) +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y="Duration of infection (days)") +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  theme_light()


ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=2, nrow=2, labels = c("A", "B", "C", "D"))


# SENSITIVITY ANALYSIS ----------------------------------------------------

### H1N1
## ferret 5775; remove "rebound" positive test at 11 dpe
H1N1_RC_ferrets[69,4] <- 0.5

## duration of infection - H1N1

alt.H1N1.infx.lengths <- data.frame()

for (ferret in H1N1_recipient_names){
  tmp.data <- H1N1_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## assume that viral growth begins immediately following the last negative test
  if (which.max(tmp.data$nw_titer > LOD) > 1){
    tmp.pos.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    tmp.pos.time <- 0
  }
  ## find the time infection resolves
  tmp.post.infection <- tmp.data %>%
    filter(dpe > tmp.pos.time)
  if (tmp.post.infection[[length(tmp.post.infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    tmp.end.time <- tmp.post.infection[[length(tmp.post.infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    tmp.end.time <- tmp.post.infection[[which.max(tmp.post.infection$nw_titer <= LOD), "dpe"]]
  }
  alt.H1N1.infx.lengths <- rbind(alt.H1N1.infx.lengths, cbind(tmp.data[1, 1:2], "duration"=tmp.end.time - tmp.pos.time))
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## linear regression for duration of infection
alt.H1N1.infx.lengths$numeric_dose <- c(6, 6, 6, 6, 4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0)
alt_H1N1_duration_regression <- lm(duration ~ numeric_dose, alt.H1N1.infx.lengths)
## no change

H1N1.infx.lengths$st <- "wt"
alt.H1N1.infx.lengths$st <- "m"
st.H1N1.infx.lengths <- rbind(H1N1.infx.lengths, alt.H1N1.infx.lengths)

alt_panel_b <- ggplot(st.H1N1.infx.lengths, aes(x=donor_dose, y=duration, shape=st)) +
  geom_beeswarm(cex=2,size=2, fill=H1N1_color, color=H1N1_color) +
  ## add regression line
  geom_abline(slope = coef(H1N1_duration_regression)[[2]], 
              intercept = coef(H1N1_duration_regression)[[1]], 
              color="black", linewidth=1) +
  geom_abline(slope = coef(alt_H1N1_duration_regression)[[2]], 
              intercept = coef(alt_H1N1_duration_regression)[[1]], 
              color="black", linewidth=1, linetype=2) +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y="Duration of infection (days)") +
  scale_shape_manual(values=c(1,19)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  guides(shape="none") +
  theme_light()

### H3N2

## ferret 626; remove initial positive test followed by a negative test
H3N2_RC_ferrets[49,4] <- 0.5
## ferret 2124; remove initial positive test followed by a negative test
H3N2_RC_ferrets[2,4] <- 0.5

## duration of infection - H3N2

alt.H3N2.infx.lengths <- data.frame()

for (ferret in H3N2_recipient_names){
  tmp.data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## assume that viral growth begins immediately following the last negative test
  if (which.max(tmp.data$nw_titer > LOD) > 1){
    tmp.pos.time <- tmp.data[[which.max(tmp.data$nw_titer > LOD)-1, "dpe"]]
  } else { ## if the first positive test is the first test, assume that
    ## the contact animal was infected immediately following exposure
    tmp.pos.time <- 0
  }
  ## find the time infection resolves
  tmp.post.infection <- tmp.data %>%
    filter(dpe > tmp.pos.time)
  if (tmp.post.infection[[length(tmp.post.infection$dpe), "nw_titer"]] > LOD){
    ## if the animal doesn't resolve infection by the last measured timepoint, 
    ## assume that infection would have resolved at the next test date 
    tmp.end.time <- tmp.post.infection[[length(tmp.post.infection$dpe), "dpe"]] + 2
  } else {
    ## if the animal does resolve, just find the first time titers reach LOD again
    tmp.end.time <- tmp.post.infection[[which.max(tmp.post.infection$nw_titer <= LOD), "dpe"]]
  }
  alt.H3N2.infx.lengths <- rbind(alt.H3N2.infx.lengths, cbind(tmp.data[1, 1:2], "duration"=tmp.end.time - tmp.pos.time))
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## linear regression for duration of infection
alt.H3N2.infx.lengths$numeric_dose <- c(6, 6, 6, 4, 4, 4, 4, 3, 2, 2, 1)
alt_H3N2_duration_regression <- lm(duration ~ numeric_dose, alt.H3N2.infx.lengths)
## no change

H3N2.infx.lengths$st <- "wt"
alt.H3N2.infx.lengths$st <- "m"
st.H3N2.infx.lengths <- rbind(H3N2.infx.lengths, alt.H3N2.infx.lengths)

alt_panel_d <- ggplot(st.H3N2.infx.lengths, aes(x=donor_dose, y=duration, shape=st)) +
  geom_beeswarm(cex=2, size=2, fill=H3N2_color, color=H3N2_color) +
  ## add regression line
  geom_abline(slope = coef(H3N2_duration_regression)[[2]], 
              intercept = coef(H3N2_duration_regression)[[1]], 
              color="black", linewidth=1) +
  geom_abline(slope = coef(alt_H3N2_duration_regression)[[2]], 
              intercept = coef(alt_H3N2_duration_regression)[[1]], 
              color="black", linewidth=1, linetype=2) +
  scale_shape_manual(values=c(1,19)) +
  labs(x=expression(paste("Index inoculum dose (", TCID[50],"/mL", ")")), y="Duration of infection (days)") +
  guides(shape = "none") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(limits = c("10^0", "10^1", "10^2", "10^3", "10^4", "10^6"), 
                   labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  theme_light()

ggarrange(panel_a, alt_panel_b, panel_c, alt_panel_d, ncol=2, nrow=2, labels = c("A", "B", "C", "D"))
