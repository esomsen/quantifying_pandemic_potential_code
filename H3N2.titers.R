library(tidyverse)
library(ggpubr)
library(khroma)

H3N2_color <- color("muted")(2)
H3N2_color <- H3N2_color[[2]]

ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "dpi", "nw_titer", "donor_dose")

## donor and contact together

all_10.6 <- ferrets %>%
  filter(Dose == "10^6" | donor_dose == "10^6") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.6$pair_shape <- c(rep(0, nrow(all_10.6)))

for (row in 1:nrow(all_10.6)){
  id <- all_10.6[row, "Ferret_ID"]
  pair <- all_10.6[row, "DI_RC_Pair"]
  if (id == "1408" | pair == "1408"){
    all_10.6$pair_shape[row] <- "a"
  } else if (id == "2124" | pair == "2124"){
    all_10.6$pair_shape[row] <- "b"
  } else if (id == "2119" | pair == "2119"){
    all_10.6$pair_shape[row] <- "c"
  } else if (id == "2104" | pair == "2104"){
    all_10.6$pair_shape[row] <- "d"
  }
}

all_10.6$shape_combo <- paste(all_10.6$LOD_shape, "-", all_10.6$pair_shape)

all_10.4 <- ferrets %>%
  filter(Dose == "10^4" | donor_dose == "10^4") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.4$pair_shape <- c(rep(0, nrow(all_10.4)))

for (row in 1:nrow(all_10.4)){
  id <- all_10.4[row, "Ferret_ID"]
  pair <- all_10.4[row, "DI_RC_Pair"]
  if (id == "F6322" | pair == "F6322"){
    all_10.4$pair_shape[row] <- "a"
  } else if (id == "F6324" | pair == "F6324"){
    all_10.4$pair_shape[row] <- "b"
  } else if (id == "F6334" | pair == "F6334"){
    all_10.4$pair_shape[row] <- "c"
  } else if (id == "F6336" | pair == "F6336"){
    all_10.4$pair_shape[row] <- "d"
  }
}

## need to remove missing titer point for now
all_10.4 <- all_10.4[-20,]

all_10.4$shape_combo <- paste(all_10.4$LOD_shape, "-", all_10.4$pair_shape)

all_10.3 <- ferrets %>%
  filter(Dose == "10^3" | donor_dose == "10^3") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.3$pair_shape <- c(rep(0, nrow(all_10.3)))

for (row in 1:nrow(all_10.3)){
  id <- all_10.3[row, "Ferret_ID"]
  pair <- all_10.3[row, "DI_RC_Pair"]
  if (id == "7439" | pair == "7439"){
    all_10.3$pair_shape[row] <- "a"
  } else if (id == "7441" | pair == "7441"){
    all_10.3$pair_shape[row] <- "b"
  } else if (id == "7425" | pair == "7425"){
    all_10.3$pair_shape[row] <- "c"
  } else if (id == "7429" | pair == "7429"){
    all_10.3$pair_shape[row] <- "d"
  }
}

all_10.3$shape_combo <- paste(all_10.3$LOD_shape, "-", all_10.3$pair_shape)

all_10.2 <- ferrets %>%
  filter(Dose == "10^2" | donor_dose == "10^2") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.2$pair_shape <- c(rep(0, nrow(all_10.2)))

for (row in 1:nrow(all_10.2)){
  id <- all_10.2[row, "Ferret_ID"]
  pair <- all_10.2[row, "DI_RC_Pair"]
  if (id == "F6314" | pair == "F6314"){
    all_10.2$pair_shape[row] <- "a"
  } else if (id == "626" | pair == "626"){
    all_10.2$pair_shape[row] <- "b"
  } else if (id == "F6328" | pair == "F6328"){
    all_10.2$pair_shape[row] <- "c"
  } else if (id == "F6330" | pair == "F6330"){
    all_10.2$pair_shape[row] <- "d"
  }
}

all_10.2$shape_combo <- paste(all_10.2$LOD_shape, "-", all_10.2$pair_shape)

all_10.1 <- ferrets %>%
  filter(Dose == "10^1" | donor_dose == "10^1") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.1$pair_shape <- c(rep(0, nrow(all_10.1)))

for (row in 1:nrow(all_10.1)){
  id <- all_10.1[row, "Ferret_ID"]
  pair <- all_10.1[row, "DI_RC_Pair"]
  if (id == "7430" | pair == "7430"){
    all_10.1$pair_shape[row] <- "a"
  } else if (id == "7434" | pair == "7434"){
    all_10.1$pair_shape[row] <- "b"
  } else if (id == "7426" | pair == "7426"){
    all_10.1$pair_shape[row] <- "c"
  } else if (id == "7" | pair == "7"){
    all_10.1$pair_shape[row] <- "d"
  }
}

all_10.1$shape_combo <- paste(all_10.1$LOD_shape, "-", all_10.1$pair_shape)

all_10.0 <- ferrets %>%
  filter(Dose == "10^0" | donor_dose == "10^0") %>%
  mutate(dpch = dpi - 1) %>%
  mutate(LOD_shape = ifelse(nw_titer <= 0.5, "below", "above")) 

all_10.0$pair_shape <- c(rep(0, nrow(all_10.1)))

for (row in 1:nrow(all_10.1)){
  id <- all_10.0[row, "Ferret_ID"]
  pair <- all_10.0[row, "DI_RC_Pair"]
  if (id == "7442" | pair == "7442"){
    all_10.0$pair_shape[row] <- "a"
  } else if (id == "7443" | pair == "7443"){
    all_10.0$pair_shape[row] <- "b"
  } else if (id == "7422" | pair == "7422"){
    all_10.0$pair_shape[row] <- "c"
  } else if (id == "7424" | pair == "7424"){
    all_10.0$pair_shape[row] <- "d"
  }
}

all_10.0$shape_combo <- paste(all_10.0$LOD_shape, "-", all_10.0$pair_shape)

p_10.6 <- ggplot(all_10.6, aes(x=dpch, y=nw_titer, color=DI_RC, group=pair_shape)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 75%", x=11, y=5, color="black") +
  labs(title = expression(paste("F; ", 10^{6})), x = NULL, y=NULL, color=NULL)

p_10.4 <- ggplot(all_10.4, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 100%", x=11, y=5, color="black") +
  labs(title = expression(paste("E; ", 10^{4})), x = NULL, y = NULL, color=NULL)

p_10.3 <- ggplot(all_10.3, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 25%", x=11, y=5, color="black") +
  labs(title = expression(paste("D; ", 10^{3})), x = NULL, y = NULL, color=NULL)

p_10.2 <- ggplot(all_10.2, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 50%", x=11, y=5, color="black") +
  labs(title = expression(paste("C; ", 10^{2})), x = NULL, y = NULL, color=NULL)

p_10.1 <- ggplot(all_10.1, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 25%", x=11, y=5, color="black") +
  labs(title = expression(paste("B; ", 10^{1})), x = NULL, y = NULL, color=NULL)

p_10.0 <- ggplot(all_10.0, aes(x=dpch, y=nw_titer, color=DI_RC)) +
  geom_point(aes(shape=shape_combo), size=2) +
  geom_line(aes(group=Ferret_ID), linewidth=1) +
  scale_color_manual(labels = c("Index", "Contact"), values = c("black", H3N2_color)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 0, 1, 2, 5)) +
  scale_x_continuous(limits=c(0, 13), breaks = seq(0, 13, 2)) +
  scale_y_continuous(limits=c(0, 7), breaks = seq(0, 7, 2)) +
  guides(shape = "none") +
  theme(legend.position = "top") +
  theme_light() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_text(label="TE = 0%", x=11, y=5, color="black") +
  labs(title = expression(paste("A; ", 10^{0})), x = NULL, y = NULL, color=NULL)

p_all <- ggarrange(p_10.0, p_10.1, p_10.2, p_10.3, p_10.4, p_10.6, 
                   ncol = 2, 
                   nrow = 3, 
                   common.legend = T, 
                   legend = "none")

p_all <- annotate_figure(p_all, left = text_grob(expression(paste("Viral titer (", log[10], TCID[50], ")")), rot = 90), bottom = "Days post exposure")

ggarrange(p_all, H3N2_kinetcs_plot, common.legend = T)
