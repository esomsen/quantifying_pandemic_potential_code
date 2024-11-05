library(khroma)
library(tidyverse)
library(ggpubr)
library(DescTools)

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

## result stays significant

panel_a <- ggplot(H3N2.init.titers, aes(x=numeric_dose, y=nw_titer)) +
  geom_point(size=2, position=position_jitter(width=0.2, height=0), fill="black", color="black") +
  ## add regression line
  geom_abline(slope = coef(H3N2.init.regression)[[2]], 
              intercept = coef(H3N2.init.regression)[[1]], 
              color="black", linewidth=1) +
  guides(color = "none") +
  ## signif
  annotate("text", x=3:3, y=7:5.5, label=c("p < 0.05", "R ^ 2 == 0.392"), parse=T, size=8, color="black") +
  labs(title="A", x=expression(paste("Index dose (", TCID[50], "/mL)")), y=expression(paste("Index initial titer (",TCID[50], "/mL)"))) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), limits = c(0, 8), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_x_continuous(breaks=seq(0, 6, 1), limits=c(-0.2, 6.2), labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2)

## linear regression for index AUC

H3N2.donor.AUC$numeric_dose <- as.numeric(substr(H3N2.donor.AUC$dose, 4, 4))
H3N2.donor.AUC.regression <- lm(AUC ~ numeric_dose, H3N2.donor.AUC)

## result becomes INSIGNIFICANT

panel_b <- ggplot(H3N2.donor.AUC, aes(x=numeric_dose, y=AUC, color=type, fill=type)) +
  geom_point(size=2, position=position_jitter(width=0.2, height=0), fill="black", color="black") +
  ## add regression line
  geom_abline(slope = coef(H3N2.donor.AUC.regression)[[2]], 
              intercept = coef(H3N2.donor.AUC.regression)[[1]], 
              color="black", linewidth=1) +
  labs(title="B", x=expression(paste("Index dose (",TCID[50], "/mL)")), y="Index AUC") +
  guides(color="none") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 4)) +
  scale_x_continuous(breaks=seq(0, 6, 1), limits=c(-0.2, 6.2), labels=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  theme_light()

## infection outcome by AUC

H3N2.AUC.infx <- merge(x=H3N2.donor.AUC[,c(1, 6)], y=H3N2_ferrets[,c(1,3, 4)], by="Ferret_ID", sort=F, all.x=T) %>%
  unique()
H3N2.AUC.infx$infx.outcome = ifelse(H3N2.AUC.infx$DI_RC_Pair %in% H3N2_recipient_names, 1, 0)

H3N2.logit <- glm(infx.outcome ~ AUC, data=H3N2.AUC.infx, family="binomial")

H3N2.AUC.infx$logit <- predict(H3N2.logit, H3N2.AUC.infx[,c(2,5)], type="response")
#save(H3N2.AUC.infx, file="H3_10_6/H3N2_no_10.6_empirical_logit.RData")

## result stays insignificant

panel_c <- ggplot(H3N2.AUC.infx, aes(x=AUC, y=infx.outcome, shape=dose)) +
  geom_point(size=2, fill=H3N2_color, color=H3N2_color) +
  scale_shape_manual(values=c(3, 16, 4, 17, 18)) +
  geom_line(aes(x=AUC, y=logit), color=H3N2_color, linewidth=1) +
  labs(title="C", x="Index AUC", y="Contact infection outcome", shape="Index dose") +
  scale_y_continuous(limits=c(0,1), breaks=c(0, 1)) +
  xlim(0, 30) +
  theme_light()

## time of initial positive test in contacts vs initial viral titer

H3N2.time.positive <- data.frame()

for (ferret in H3N2_recipient_names){
  ferret_data <- H3N2_RC_ferrets %>%
    filter(Ferret_ID == ferret)
  ## find first positive test
  first_positive_time <- ferret_data[[which.max(ferret_data$nw_titer > LOD), "dpe"]]
  H3N2.time.positive <- rbind(H3N2.time.positive, c(ferret_data[1,1:2], "pos_time" = first_positive_time))
}

H3N2.time.positive <- merge(x=H3N2.time.positive, y=H3N2_ferrets[,c(1,3)], by.x = "Ferret_ID", by.y="DI_RC_Pair", sort=F, all.x=T) %>%
  unique()
H3N2.time.positive <- merge(x=H3N2.time.positive, y=H3N2.init.titers[,c(1,4)], by.x="Ferret_ID.y", by.y="Ferret_ID", sort=F, allx.=T) %>%
  unique()
#H3N2.time.positive$donor_dose <- as.factor(substr(H3N2.time.positive$donor_dose, 4, 4))
names(H3N2.time.positive) <- c("index", "contact", "donor_dose", "pos_time", "initial_titer")

## linear regression for time to first positive test
H3N2.time.positive.regression <- lm(pos_time ~ initial_titer, H3N2.time.positive)

## result stays insignificant

panel_d <- ggplot(H3N2.time.positive, aes(x=initial_titer, y=pos_time, shape=donor_dose)) +
  geom_point(size=2, position=position_jitter(width=0.2, height=0), fill=H3N2_color, color=H3N2_color) +
  scale_shape_manual(values=c(3, 16, 4, 17, 18), labels=c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^6))) +
  ## add regression line
  geom_abline(slope = coef(H3N2.time.positive.regression)[[2]], 
              intercept = coef(H3N2.time.positive.regression)[[1]], 
              color=H3N2_color, linewidth=1) +
  labs(title="D", x=expression(paste("Initial titer of index (",TCID[50], "/mL)")), y="Time of first positive test (days)", shape="Index dose") +
  guides(color = "none") +
  scale_x_continuous(limits = c(0, 8), breaks=seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  theme_light()

ggarrange(panel_a, panel_b, panel_c, panel_d, ncol=4, align="h", common.legend = T, legend="right")
