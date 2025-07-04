library(tidyverse)
library(DescTools)

LOD <- 1

interpolation <- function(row1, row2, data, interval=0.1){
  fxn.times <- seq(data[[row1,"time"]], data[[row2,"time"]], interval)
  fxn.preds <- seq(data[[row1,"nw_titer"]], data[[row2,"nw_titer"]], length.out=length(fxn.times))
  fxn.df <- data.frame(time = fxn.times, titer = fxn.preds)
  return (fxn.df)
  rm(list=ls(pattern="^fxn"))
}


## import data
H1N1_ferrets <- read_csv("H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H1N1_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "time", "nw_titer", "donor_dose")

H1N1_DI_ferrets <- H1N1_ferrets %>%
  ## keep only donor ferrets
  filter(DI_RC == "DI") %>%
  dplyr::select(Ferret_ID, Dose, time, nw_titer) %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  ## create column for days post exposure
  mutate(dpe = time - 1)
H1N1_donor_names <- unique(H1N1_DI_ferrets$Ferret_ID)

H3N2_ferrets <- read_csv("H3N2_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(H3N2_ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "time", "nw_titer", "donor_dose")
## removing tests at 13 and 14 dpi because we assume that there is no impact on transmission for tests at or below LOD
H3N2_ferrets <- H3N2_ferrets %>%
  filter(!time %in% c(13, 14))

H3N2_DI_ferrets <- H3N2_ferrets %>%
  ## keep only donor ferrets
  filter(DI_RC == "DI") %>%
  dplyr::select(Ferret_ID, Dose, time, nw_titer) %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  ## create column for days post exposure
  mutate(dpe = time - 1)

## removing F6335 because of missing data point
## removing 7423 because this animal doesn't get infected
## remove all donors who never became infected, plus their pair (the 10^0 animals)
H3N2_DI_ferrets <- H3N2_DI_ferrets %>%
  filter(Ferret_ID != "F6335") %>%
  filter(Ferret_ID != "7423") %>%
  filter(Dose != "10^0")

H3N2_donor_names <- unique(H3N2_DI_ferrets$Ferret_ID)

## create list to store interpolated titers
H1N1.ferret.preds <- vector("list", length(H1N1_donor_names))
names(H1N1.ferret.preds) <- H1N1_donor_names

for (ferret in H1N1_donor_names){
  ## extract data for each donor
  tmp.data <- H1N1_DI_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    select(c(time, nw_titer))
  ## create temporary df to store interpolations
  tmp.df <- data.frame()
  ## loop through all data until last timepoint
  for (t in 1:(length(tmp.data$time)-1)){
    ## interpolation
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data))
  }
  ## force numeric, remove duplicate rows
  tmp.df <- tmp.df %>%
    distinct() %>%
    mutate(time = as.numeric(time))
  ## store interpolated titers in list
  H1N1.ferret.preds[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## create list to store interpolated titers
H3N2.ferret.preds <- vector("list", length(H3N2_donor_names))
names(H3N2.ferret.preds) <- H3N2_donor_names

for (ferret in H3N2_donor_names){
  ## extract data for each donor
  tmp.data <- H3N2_DI_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    select(c(time, nw_titer))
  ## create temporary df to store interpolations
  tmp.df <- data.frame()
  ## loop through all data until last timepoint
  for (t in 1:(length(tmp.data$time)-1)){
    ## interpolation
    tmp.df <- rbind(tmp.df, interpolation(t, t+1, tmp.data))
  }
  ## force numeric, remove duplicate rows
  tmp.df <- tmp.df %>%
    distinct() %>%
    mutate(time = as.numeric(time))
  ## store interpolated titers in list
  H3N2.ferret.preds[[ferret]] <- tmp.df
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
}

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## range of s values to try
s_vals <- seq(0.001, 0.2, 0.001)
## log prob tracker
H1N1.log.probs <- c()

## for H1N1
for (s in s_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H1N1_donor_names)))
  names(tmp.probs) <- H1N1_donor_names
  
  for (ferret in H1N1_donor_names){
    tmp.contact <- H1N1_ferrets %>%
      filter(DI_RC == "RC") %>%
      filter(DI_RC_Pair == ferret)
    ## if contact is infected, use this to calc prob of s
    if (max(tmp.contact$nw_titer) > LOD){
      ## find time of first positive test and last negative test in paired contact
      tmp.pos.time <- tmp.contact[[which.max(tmp.contact$nw_titer > LOD),"time"]]
      ## if the first positive test is at 2dpi, call the last negative time at 0
      ## by definition, the integral at time 0 is 0
      if (tmp.pos.time == 2){
        tmp.neg.time <- 0
        tmp.neg.integral <- 0
      } else {
        tmp.neg.time <- tmp.contact[[which.max(tmp.contact$nw_titer > LOD)-1, "time"]]
        ## find integral value at time of last negative test
        tmp.neg.titers <- H1N1.ferret.preds[[ferret]] %>%
          filter(time <= tmp.neg.time)
        tmp.yvals <- tmp.neg.titers$titer
        ## when titer is at or below LOD, AUC is 0
        tmp.below_LOD <- which(tmp.neg.titers$titer < LOD)
        if (length(tmp.below_LOD) >= 1){
          tmp.yvals[tmp.below_LOD] <- 0
        }
        tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.yvals*s, method="trapezoid")
      }
      ## find integral from last neg to first positive test
      tmp.df <- H1N1.ferret.preds[[ferret]] %>%
        filter(time >= tmp.neg.time) %>%
        filter(time <= tmp.pos.time)
      tmp.yvals <- tmp.df$titer
      ## when titer is at or below LOD, AUC is 0
      tmp.below_LOD <- which(tmp.df$titer < LOD)
      if (length(tmp.below_LOD) >= 1){
        tmp.yvals[tmp.below_LOD] <- 0
      }
      tmp.AUC <- AUC(x=tmp.df$time, y=tmp.yvals*s, method="trapezoid")
      ## find prob of getting infected between last neg and first pos test
      tmp.interval.prob <- 1 - exp(-tmp.AUC)
      ## overall prob is 
      ## prob of not getting infected between 1dpi and last negative test 
      ## TIMES interval prob
      tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
    } else { ## if contact is uninfected, use this to calc prob of s
      tmp.yvals <- H1N1.ferret.preds[[ferret]]$titer
      ## when titer is at or below LOD, AUC is 0
      tmp.below_LOD <- which(H1N1.ferret.preds[[ferret]]$titer < LOD)
      if (length(tmp.below_LOD) >= 1){
        tmp.yvals[tmp.below_LOD] <- 0
      }
      tmp.AUC <- AUC(x=H1N1.ferret.preds[[ferret]]$time, y=tmp.yvals*s, method="trapezoid")
      tmp.prob <- exp(-tmp.AUC)
    }
    tmp.probs[ferret] <- tmp.prob
  }
  tmp.log.prs <- log(tmp.probs)
  ## replace any Inf log probs with log of a really small number
  tmp.log.prs <- replace(tmp.log.prs, tmp.log.prs==0, -744.4401)
  tmp.sum.prs <- sum(tmp.log.prs)
  H1N1.log.probs <- append(H1N1.log.probs, tmp.sum.prs)
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", s, "done"))
}

## log prob tracker
H3N2.log.probs <- c()

## for H3N2
for (s in s_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H3N2_donor_names)))
  names(tmp.probs) <- H3N2_donor_names
  
  for (ferret in H3N2_donor_names){
    tmp.contact <- H3N2_ferrets %>%
      filter(DI_RC == "RC") %>%
      filter(DI_RC_Pair == ferret)
    ## if contact is infected, use this to calc prob of s
    if (max(tmp.contact$nw_titer) > LOD){
      ## find time of first positive test and last negative test in paired contact
      tmp.pos.time <- tmp.contact[[which.max(tmp.contact$nw_titer > LOD),"time"]]
      ## if the first positive test is at 2dpi, call the last negative time at 0
      ## by definition, the integral at time 0 is 0
      if (tmp.pos.time == 2){
        tmp.neg.time <- 0
        tmp.neg.integral <- 0
      } else {
        tmp.neg.time <- tmp.contact[[which.max(tmp.contact$nw_titer > LOD)-1, "time"]]
        ## find integral value at time of last negative test
        tmp.neg.titers <- H3N2.ferret.preds[[ferret]] %>%
          filter(time <= tmp.neg.time)
        tmp.yvals <- tmp.neg.titers$titer
        ## when titer is at or below LOD, AUC is 0
        tmp.below_LOD <- which(tmp.neg.titers$titer < LOD)
        if (length(tmp.below_LOD) >= 1){
          tmp.yvals[tmp.below_LOD] <- 0
        }
        tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.yvals*s, method="trapezoid")
      }
      ## find integral from last neg to first positive test
      tmp.df <- H3N2.ferret.preds[[ferret]] %>%
        filter(time >= tmp.neg.time) %>%
        filter(time <= tmp.pos.time)
      tmp.yvals <- tmp.df$titer
      ## when titer is at or below LOD, AUC is 0
      tmp.below_LOD <- which(tmp.df$titer < LOD)
      if (length(tmp.below_LOD) >= 1){
        tmp.yvals[tmp.below_LOD] <- 0
      }
      tmp.AUC <- AUC(x=tmp.df$time, y=tmp.yvals*s, method="trapezoid")
      ## find prob of getting infected between last neg and first pos test
      tmp.interval.prob <- 1 - exp(-tmp.AUC)
      ## overall prob is 
      ## prob of not getting infected between 1dpi and last negative test 
      ## TIMES interval prob
      tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
    } else { ## if contact is uninfected, use this to calc prob of s
      tmp.yvals <- H3N2.ferret.preds[[ferret]]$titer
      ## when titer is at or below LOD, AUC is 0
      tmp.below_LOD <- which(H3N2.ferret.preds[[ferret]]$titer < LOD)
      if (length(tmp.below_LOD) >= 1){
        tmp.yvals[tmp.below_LOD] <- 0
      }
      tmp.AUC <- AUC(x=H3N2.ferret.preds[[ferret]]$time, y=tmp.yvals*s, method="trapezoid")
      tmp.prob <- exp(-tmp.AUC)
    }
    tmp.probs[ferret] <- tmp.prob
  }
  tmp.log.prs <- log(tmp.probs)
  ## replace any Inf log probs with log of a really small number
  tmp.log.prs <- replace(tmp.log.prs, tmp.log.prs==0, -744.4401)
  tmp.sum.prs <- sum(tmp.log.prs)
  H3N2.log.probs <- append(H3N2.log.probs, tmp.sum.prs)
  ## tidy environment
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", s, "done"))
}

prob.trace <- data.frame(s = s_vals,
                         H1N1.pr = H1N1.log.probs,
                         H3N2.pr = H3N2.log.probs)

H1N1.MLE <- prob.trace[which.max(H1N1.log.probs),1]
H1N1.MLE.logL <- prob.trace[which.max(H1N1.log.probs),2]
H3N2.MLE <- prob.trace[which.max(H3N2.log.probs),1]
H3N2.MLE.logL <- prob.trace[which.max(H3N2.log.probs),3] 

H1N1.CIs <- prob.trace[near(prob.trace[which.max(H1N1.log.probs),2]-1.92, prob.trace$H1N1.pr, tol=0.11),1:2]
H1N1.CIs <- H1N1.CIs[c(2,5),1]

H3N2.CIs <- prob.trace[near(prob.trace[which.max(H3N2.log.probs),3]-1.92, prob.trace$H3N2.pr, tol=0.11),c(1,3)]
H3N2.CIs <- H3N2.CIs[2:3,1]

save(prob.trace, file="prob.trace.Rdata")
save(H1N1.MLE.logL, file="H1N1.MLE.logL.Rdata")
save(H3N2.MLE.logL, file="H3N2.MLE.logL.Rdata")

# plots -------------------------------------------------------------------

library(khroma)
library(ggpubr)

plot_colors <- color("muted")(2)
prob.trace <- prob.trace %>%
  pivot_longer(-s, names_to="virus", values_to = "likelihood")

panel_a <- ggplot(prob.trace, aes(x=s, y=likelihood, color=virus)) +
  geom_line(linewidth=1.5) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  labs(x="Transmissibility, s", y="Log likelihood", color="Virus") +
  geom_vline(xintercept=H1N1.MLE, color=plot_colors[1], linewidth=1, linetype = 2) +
  geom_vline(xintercept=H3N2.MLE, color=plot_colors[2], linewidth=1, linetype = 2) +
  annotate("rect", xmin=H1N1.CIs[1], xmax=H1N1.CIs[2], ymin=-90, ymax=-15, alpha=0.2, color=plot_colors[1], fill=plot_colors[1]) +
  annotate("rect", xmin=H3N2.CIs[1], xmax=H3N2.CIs[2], ymin=-90, ymax=-15, alpha=0.2, color=plot_colors[2], fill=plot_colors[2]) +
  scale_y_continuous(limits=c(-90, -15), expand=c(0, 0)) +
  theme(legend.position = "top") +
  theme_light()

calculate_pr_constant <- function(s, VL){
  lambda <- s * VL
  ## integrate from 0 to 1 for one-hour exposure
  integral <- AUC(x=c(0,1), y=c(lambda, lambda), method = "trapezoid")
  prob <- 1 - exp(-integral)
  return (prob)
}

VL_probs_H1N1 <- data.frame(log_VL = seq(0, 10, 0.1), 
                            prob = rep(0, 101))

for (row in 1:nrow(VL_probs_H1N1)){
  slice <- VL_probs_H1N1[row,]
  if (slice$log_VL > LOD){
    prob <- 1 - exp(-slice$log_VL * H1N1.MLE)
  } else {
    prob <- 0
  }
  VL_probs_H1N1$prob[row] <- prob
}

VL_probs_H3N2 <- data.frame(log_VL = seq(0, 10, 0.1), 
                            prob = rep(0, 101))
for (row in 1:nrow(VL_probs_H3N2)){
  slice <- VL_probs_H3N2[row,]
  if (slice$log_VL > LOD){
    prob <- 1 - exp(-slice$log_VL * H3N2.MLE)
  } else {
    prob <- 0
  }
  VL_probs_H3N2$prob[row] <- prob
}

combo_VL_pr <- merge(VL_probs_H1N1, VL_probs_H3N2, by="log_VL", suffixes = c(".H1N1", ".H3N2"))
combo_VL_pr <- combo_VL_pr %>%
  pivot_longer(-log_VL, names_to="virus", names_prefix = "prob.", values_to = "pr")

## now, calculate CIs for ribbon
H1N1_ribbon <- data.frame(log_VL = seq(0, 10, 0.1), 
                          prob_upper = rep(0, 101), 
                          prob_lower = rep(0, 101))
for (row in 1:nrow(H1N1_ribbon)){
  slice <- H1N1_ribbon[row,]
  if (slice$log_VL > LOD){
    prob_upper <- 1 - exp(-slice$log_VL * H1N1.CIs[2])
    prob_lower <- 1 - exp(-slice$log_VL * H1N1.CIs[1])
  } else {
    prob_upper <- 0
    prob_lower <- 0
  }
  H1N1_ribbon$prob_upper[row] <- prob_upper
  H1N1_ribbon$prob_lower[row] <- prob_lower
}

H3N2_ribbon <- data.frame(log_VL = seq(0, 10, 0.1), 
                          prob_upper = rep(0, 101), 
                          prob_lower = rep(0, 101))
for (row in 1:nrow(H3N2_ribbon)){
  slice <- H3N2_ribbon[row,]
  if (slice$log_VL > LOD){
    prob_upper <- 1 - exp(-slice$log_VL * H3N2.CIs[2])
    prob_lower <- 1 - exp(-slice$log_VL * H3N2.CIs[1])
  } else {
    prob_upper <- 0
    prob_lower <- 0
  }
  H3N2_ribbon$prob_upper[row] <- prob_upper
  H3N2_ribbon$prob_lower[row] <- prob_lower
}

panel_b <- ggplot(combo_VL_pr, aes(x=log_VL, y=pr, color=virus)) +
  geom_line(linewidth=1.5) +
  geom_ribbon(data=H1N1_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[1]], alpha = 0.2, inherit.aes = F) +
  geom_ribbon(data=H3N2_ribbon, aes(x=log_VL, ymin=prob_lower, ymax=prob_upper), fill=plot_colors[[2]], alpha = 0.2, inherit.aes = F) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]])) + 
  scale_x_continuous(limits=c(0, 7), breaks = seq(0, 7, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6))) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission", color="Virus") +
  theme_light() +
  guides(color="none")

ggarrange(panel_a, panel_b, ncol=2, labels = c("A", "B"), common.legend = T, legend = "bottom")
