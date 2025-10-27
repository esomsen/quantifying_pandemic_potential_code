library(tidyverse)
library(ggpubr)
library(khroma)

H1N1_color <- color("muted")(1)
H1N1_color <- H1N1_color[[1]]

LOD <- 1

ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "dpi", "nw_titer", "donor_dose")

## donor and contact together

all_10.6 <- ferrets %>%
  filter(Dose == "10^6" | donor_dose == "10^6") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer < LOD, "below", "above")) 

all_10.6$pair_shape <- c(rep(0, nrow(all_10.6)))

for (row in 1:nrow(all_10.6)){
  id <- all_10.6[row, "Ferret_ID"]
  pair <- all_10.6[row, "DI_RC_Pair"]
  if (id == "5109" | pair == "5109"){
    all_10.6$pair_shape[row] <- "a"
  } else if (id == "5343" | pair == "5343"){
    all_10.6$pair_shape[row] <- "b"
  } else if (id == "5429" | pair == "5429"){
    all_10.6$pair_shape[row] <- "c"
  } else if (id == "5439" | pair == "5439"){
    all_10.6$pair_shape[row] <- "d"
  }
}

all_10.6$shape_combo <- paste(all_10.6$LOD_shape, "-", all_10.6$pair_shape)

all_10.4 <- ferrets %>%
  filter(Dose == "10^4" | donor_dose == "10^4") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer < LOD, "below", "above")) 

all_10.4$pair_shape <- c(rep(0, nrow(all_10.4)))

for (row in 1:nrow(all_10.4)){
  id <- all_10.4[row, "Ferret_ID"]
  pair <- all_10.4[row, "DI_RC_Pair"]
  if (id == "5445" | pair == "5445"){
    all_10.4$pair_shape[row] <- "a"
  } else if (id == "5448" | pair == "5448"){
    all_10.4$pair_shape[row] <- "b"
  } else if (id == "5433" | pair == "5433"){
    all_10.4$pair_shape[row] <- "c"
  } else if (id == "5441" | pair == "5441"){
    all_10.4$pair_shape[row] <- "d"
  }
}

all_10.4$shape_combo <- paste(all_10.4$LOD_shape, "-", all_10.4$pair_shape)

all_10.2 <- ferrets %>%
  filter(Dose == "10^2" | donor_dose == "10^2") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer < LOD, "below", "above")) 

all_10.2$pair_shape <- c(rep(0, nrow(all_10.2)))

for (row in 1:nrow(all_10.2)){
  id <- all_10.2[row, "Ferret_ID"]
  pair <- all_10.2[row, "DI_RC_Pair"]
  if (id == "5775" | pair == "5775"){
    all_10.2$pair_shape[row] <- "a"
  } else if (id == "5777" | pair == "5777"){
    all_10.2$pair_shape[row] <- "b"
  } else if (id == "5783" | pair == "5783"){
    all_10.2$pair_shape[row] <- "c"
  } else if (id == "5787" | pair == "5787"){
    all_10.2$pair_shape[row] <- "d"
  }
}

all_10.2$shape_combo <- paste(all_10.2$LOD_shape, "-", all_10.2$pair_shape)

all_10.1 <- ferrets %>%
  filter(Dose == "10^1" | donor_dose == "10^1") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer < LOD, "below", "above")) 

all_10.1$pair_shape <- c(rep(0, nrow(all_10.1)))

for (row in 1:nrow(all_10.1)){
  id <- all_10.1[row, "Ferret_ID"]
  pair <- all_10.1[row, "DI_RC_Pair"]
  if (id == "598" | pair == "598"){
    all_10.1$pair_shape[row] <- "a"
  } else if (id == "599" | pair == "599"){
    all_10.1$pair_shape[row] <- "b"
  } else if (id == "626" | pair == "626"){
    all_10.1$pair_shape[row] <- "c"
  } else if (id == "7826" | pair == "7826"){
    all_10.1$pair_shape[row] <- "d"
  }
}

all_10.1$shape_combo <- paste(all_10.1$LOD_shape, "-", all_10.1$pair_shape)

all_10.0 <- ferrets %>%
  filter(Dose == "10^0" | donor_dose == "10^0") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer < LOD, "below", "above")) 

all_10.0$pair_shape <- c(rep(0, nrow(all_10.1)))

for (row in 1:nrow(all_10.1)){
  id <- all_10.0[row, "Ferret_ID"]
  pair <- all_10.0[row, "DI_RC_Pair"]
  if (id == "7818" | pair == "7818"){
    all_10.0$pair_shape[row] <- "a"
  } else if (id == "2684" | pair == "2684"){
    all_10.0$pair_shape[row] <- "b"
  } else if (id == "7823" | pair == "7823"){
    all_10.0$pair_shape[row] <- "c"
  } else if (id == "7824" | pair == "7824"){
    all_10.0$pair_shape[row] <- "d"
  }
}

all_10.0$shape_combo <- paste(all_10.0$LOD_shape, "-", all_10.0$pair_shape)

p_10.6 <- ggplot(all_10.6, aes(x=dpch, y=nw_titer, color=DI_RC, group=pair_shape)) +
  geom_point(aes(shape=shape_combo), size=3) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("H1N1 Index", "H1N1 Contact"), values = c("black", H1N1_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  guides(shape = "none") +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2) +
  geom_text(label="TE = 4/4", x=11, y=5, color="black") +
  labs(x = NULL, y=NULL, color=NULL)

p_10.4 <- ggplot(all_10.4, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=3) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("H1N1 Index", "H1N1 Contact"), values = c("black", H1N1_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2) +
  geom_text(label="TE = 4/4", x=11, y=5, color="black") +
  labs(x = NULL, y = NULL, color=NULL)

p_10.2 <- ggplot(all_10.2, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=3) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("H1N1 Index", "H1N1 Contact"), values = c("black", H1N1_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2) +
  geom_text(label="TE = 4/4", x=11, y=5, color="black") +
  labs(x = NULL, y = NULL, color=NULL)

p_10.1 <- ggplot(all_10.1, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=3) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("H1N1 Index", "H1N1 Contact"), values = c("black", H1N1_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2) +
  geom_text(label="TE = 4/4", x=11, y=5, color="black") +
  labs(x = NULL, y = NULL, color=NULL)

p_10.0 <- ggplot(all_10.0, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=3) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("H1N1 Index", "H1N1 Contact"), values = c("black", H1N1_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = LOD, linetype = 2) +
  geom_text(label="TE = 3/4", x=11, y=5, color="black") +
  labs(x = NULL, y = NULL, color=NULL)

## create empty plot
spacer <- ggplot() +
  theme_void() +
  geom_text(aes(0, 0, label="Not performed.")) +
  xlab(NULL)

H1N1 <- ggarrange(p_10.0, p_10.1, p_10.2, spacer, p_10.4, p_10.6, 
          ncol = 1, 
          nrow = 6,
          labels = c("A", "B", "C", "D", "E", "F"), 
          common.legend = T, 
          legend = "none",
          vjust=-.1)

H1N1 <- annotate_figure(H1N1, left = text_grob(expression(paste("Viral titer (", TCID[50],"/mL", ")")), rot = 90), bottom = "Days post exposure", top="Cal/2009")
