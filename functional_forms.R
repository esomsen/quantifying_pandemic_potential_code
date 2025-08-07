library(tidyverse)
library(DescTools)
library(khroma)
library(ggpubr)

plot_colors <- color("muted")(2)

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


# LINEAR FORM -------------------------------------------------------------

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## range of s values to try
s_vals <- seq(0.000000001, 0.0001, 0.0000001)
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
        tmp.yvals <- 10^(tmp.neg.titers$titer)
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
      tmp.yvals <- 10^(tmp.df$titer)
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
      tmp.yvals <- 10^(H1N1.ferret.preds[[ferret]]$titer)
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
  ## replace any 0 prob with a very small number
  tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
  tmp.log.prs <- log(tmp.probs)
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
        tmp.yvals <- 10^(tmp.neg.titers$titer)
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
      tmp.yvals <- 10^(tmp.df$titer)
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
      tmp.yvals <- 10^(H3N2.ferret.preds[[ferret]]$titer)
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
  ## replace any 0 prob with a very small number
  tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
  tmp.log.prs <- log(tmp.probs)
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

H1N1.CIs <- prob.trace[near(prob.trace[which.max(H1N1.log.probs),2]-1.92, prob.trace$H1N1.pr, tol=0.5),1:2]
H1N1.CIs <- H1N1.CIs[c(2,3),1]

H3N2.CIs <- prob.trace[near(prob.trace[which.max(H3N2.log.probs),3]-1.92, prob.trace$H3N2.pr, tol=0.05),c(1,3)]
H3N2.CIs <- H3N2.CIs[c(2,4),1]

## plots

VLs <- seq(0, 8, 0.1)
LOD <- 1

H1.linear.lambda <- c()
H1.linear.probs <- c()
H1.linear.lower <- c()
H1.linear.upper <- c()
H3.linear.lambda <- c()
H3.linear.probs <- c()
H3.linear.lower <- c()
H3.linear.upper <- c()

for (v in VLs){
  if (v < LOD){
    H1.linear.lambda <- append(H1.linear.lambda, 0)
    H3.linear.lambda <- append(H3.linear.lambda, 0)
    H1.linear.probs <- append(H1.linear.probs, 0)
    H1.linear.lower <- append(H1.linear.lower, 0)
    H1.linear.upper <- append(H1.linear.upper, 0)
    H3.linear.probs <- append(H3.linear.probs, 0)
    H3.linear.lower <- append(H3.linear.lower, 0)
    H3.linear.upper <- append(H3.linear.upper, 0)
  } else {
    H1.linear.lambda <- append(H1.linear.lambda, AUC(x=c(0,1), y=c((10^v)*1.601e-06, (10^v)*1.601e-06), method="trapezoid"))
    H3.linear.lambda <- append(H3.linear.lambda, AUC(x=c(0,1), y=c((10^v)*1.99e-5, (10^v)*1.99e-5), method="trapezoid"))
    H1.linear.probs <- append(H1.linear.probs, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*1.601e-06, (10^v)*1.601e-06), method="trapezoid"))))
    H1.linear.lower <- append(H1.linear.lower, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*1.001e-6, (10^v)*1.001e-6), method="trapezoid"))))
    H1.linear.upper <- append(H1.linear.upper, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*2.401e-6, (10^v)*2.401e-6), method="trapezoid"))))
    H3.linear.probs <- append(H3.linear.probs, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*1.99e-5, (10^v)*1.99e-5), method="trapezoid"))))
    H3.linear.lower <- append(H3.linear.lower, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*1e-5, (10^v)*1e-5), method="trapezoid"))))
    H3.linear.upper <- append(H3.linear.upper, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*3.49e-05, (10^v)*3.49e-05), method="trapezoid"))))
  }
}

linear.lambdas <- data.frame(VL = VLs, 
                             lambda = c(H1.linear.lambda, H3.linear.lambda), 
                             Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

linear.lambda.plot <- ggplot(linear.lambdas, aes(x=VL, y=lambda, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 10)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Force of infection") +
  theme_classic()

linear.probs <- data.frame(VL = VLs, 
                           prob = c(H1.linear.probs, H3.linear.probs), 
                           lower = c(H1.linear.lower, H3.linear.lower), 
                           upper = c(H3.linear.upper, H3.linear.upper),
                           Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

linear.prob.plot <- ggplot(linear.probs, aes(x=VL, y=prob, color=Virus)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Virus, alpha=0.3)) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  guides(color="none", fill="none", alpha="none") +
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()


# THRESHOLD ---------------------------------------------------------------

s_vals <- seq(0.001, 1, 0.001)
h_vals <- seq(LOD, 7, 0.01) ## where h is threshold point
## log prob tracker
H1N1.joint.log.probs <- matrix(ncol=length(s_vals), nrow=length(h_vals))

## nested loop
for (s in s_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H1N1_donor_names)))
  names(tmp.probs) <- H1N1_donor_names
  
  for (h in h_vals){
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
          ## when titer is below h value, FOI is 0
          tmp.yvals[which(tmp.neg.titers$titer < h)] <- 0
          ## when titer is at or above h, FOI is s
          tmp.yvals[which(tmp.neg.titers$titer >= h)] <- s
          tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.yvals, method="trapezoid")
        }
        ## find integral from last neg to first positive test
        tmp.df <- H1N1.ferret.preds[[ferret]] %>%
          filter(time >= tmp.neg.time) %>%
          filter(time <= tmp.pos.time)
        tmp.yvals <- tmp.df$titer
        ## when titer is below h value, FOI is 0
        tmp.yvals[which(tmp.df$titer < h)] <- 0
        ## when titer is at or above h, FOI is s
        tmp.yvals[which(tmp.df$titer >= h)] <- s
        tmp.AUC <- AUC(x=tmp.df$time, y=tmp.yvals, method="trapezoid")
        ## find prob of getting infected between last neg and first pos test
        tmp.interval.prob <- 1 - exp(-tmp.AUC)
        ## overall prob is 
        ## prob of not getting infected between 1dpi and last negative test 
        ## TIMES interval prob
        tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
      } else { ## if contact is uninfected, use this to calc prob of s
        tmp.yvals <- H1N1.ferret.preds[[ferret]]$titer
        ## when titer is at or below h, AUC is 0
        tmp.yvals[which(H1N1.ferret.preds[[ferret]]$titer < h)] <- 0
        ## when titer is at or above h, FOI is s
        tmp.yvals[which(H1N1.ferret.preds[[ferret]]$titer >= h)] <- s
        tmp.AUC <- AUC(x=H1N1.ferret.preds[[ferret]]$time, y=tmp.yvals, method="trapezoid")
        tmp.prob <- exp(-tmp.AUC)
      }
      tmp.probs[ferret] <- tmp.prob
    }
    ## replace any 0 prob with a very small number
    tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
    tmp.sum.prs <- sum(log(tmp.probs))
    ## for each s value, record likelihood of all h values
    H1N1.joint.log.probs[which(h_vals == h), which(s_vals == s)] <- tmp.sum.prs
  }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", s, "done"))
}  

## log prob tracker
H3N2.joint.log.probs <- matrix(ncol=length(s_vals), nrow=length(h_vals))

## nested loop
for (s in s_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H3N2_donor_names)))
  names(tmp.probs) <- H3N2_donor_names
  
  for (h in h_vals){
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
          ## when titer is below h value, FOI is 0
          tmp.yvals[which(tmp.neg.titers$titer < h)] <- 0
          ## when titer is at or above h, FOI is s
          tmp.yvals[which(tmp.neg.titers$titer >= h)] <- s
          tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.yvals, method="trapezoid")
        }
        ## find integral from last neg to first positive test
        tmp.df <- H3N2.ferret.preds[[ferret]] %>%
          filter(time >= tmp.neg.time) %>%
          filter(time <= tmp.pos.time)
        tmp.yvals <- tmp.df$titer
        ## when titer is below h value, FOI is 0
        tmp.yvals[which(tmp.df$titer < h)] <- 0
        ## when titer is at or above h, FOI is s
        tmp.yvals[which(tmp.df$titer >= h)] <- s
        tmp.AUC <- AUC(x=tmp.df$time, y=tmp.yvals, method="trapezoid")
        ## find prob of getting infected between last neg and first pos test
        tmp.interval.prob <- 1 - exp(-tmp.AUC)
        ## overall prob is 
        ## prob of not getting infected between 1dpi and last negative test 
        ## TIMES interval prob
        tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
      } else { ## if contact is uninfected, use this to calc prob of s
        tmp.yvals <- H3N2.ferret.preds[[ferret]]$titer
        ## when titer is at or below h, AUC is 0
        tmp.yvals[which(H3N2.ferret.preds[[ferret]]$titer < h)] <- 0
        ## when titer is at or above h, FOI is s
        tmp.yvals[which(H3N2.ferret.preds[[ferret]]$titer >= h)] <- s
        tmp.AUC <- AUC(x=H3N2.ferret.preds[[ferret]]$time, y=tmp.yvals, method="trapezoid")
        tmp.prob <- exp(-tmp.AUC)
      }
      tmp.probs[ferret] <- tmp.prob
    }
    ## replace any 0 prob with a very small number
    tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
    tmp.sum.prs <- sum(log(tmp.probs))
    ## for each s value, record likelihood of all h values
    H3N2.joint.log.probs[which(h_vals == h), which(s_vals == s)] <- tmp.sum.prs
  }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", s, "done"))
}  

which(H1N1.joint.log.probs == max(H1N1.joint.log.probs), arr.ind=T)
H1N1.MLE <- c("s" = s_vals[708], "h" = h_vals[284])
which(near(H1N1.joint.log.probs,H1N1.joint.log.probs[284,708]-1.92, tol=0.00008), arr.ind=T)

which(H3N2.joint.log.probs == max(H3N2.joint.log.probs), arr.ind=T)
H3N2.MLE <- c("s" = s_vals[118], "h" = h_vals[1])
H3N2.joint.log.probs[1,118]

## plots

H1.threshold.lambda <- c()
H1.threshold.probs <- c()
H1.threshold.lower <- c()
H1.threshold.upper <- c()
H3.threshold.lambda <- c()
H3.threshold.probs <- c()
H3.threshold.lower <- c()
H3.threshold.upper <- c()

for (v in VLs){
  if (v < 3.83){
    H1.threshold.lambda <- append(H1.threshold.lambda, 0)
    H1.threshold.probs <- append(H1.threshold.probs, 0)
    H1.threshold.lower <- append(H1.threshold.lower, 0)
    H1.threshold.upper <- append(H1.threshold.upper, 0)
  } else {
    H1.threshold.lambda <- append(H1.threshold.lambda, AUC(x=c(0, 1), y=c(0.708, 0.708), method="trapezoid"))
    H1.threshold.probs <- append(H1.threshold.probs, (1 - exp(-AUC(x=c(0, 1), y=c(0.708, 0.708), method="trapezoid"))))
    H1.threshold.lower <- append(H1.threshold.lower, (1 - exp(-AUC(x=c(0, 1), y=c(0.426, 0.426), method="trapezoid"))))
    H1.threshold.upper <- append(H1.threshold.upper, (1 - exp(-AUC(x=c(0, 1), y=c(1.112, 1.112), method="trapezoid"))))
  }
}

for (v in VLs){
  if (v < 1){
    H3.threshold.lambda <- append(H3.threshold.lambda, 0)
    H3.threshold.probs <- append(H3.threshold.probs, 0)
    H3.threshold.lower <- append(H3.threshold.lower, 0)
    H3.threshold.upper <- append(H3.threshold.upper, 0)
  } else {
    H3.threshold.lambda <- append(H3.threshold.lambda, AUC(x=c(0, 1), y=c(0.118, 0.118), method="trapezoid"))
    H3.threshold.probs <- append(H3.threshold.probs, (1 - exp(-AUC(x=c(0, 1), y=c(0.118, 0.118), method="trapezoid"))))
    H3.threshold.lower <- append(H3.threshold.lower, (1 - exp(-AUC(x=c(0, 1), y=c(0.06, 0.06), method="trapezoid"))))
    H3.threshold.upper <- append(H3.threshold.upper, (1 - exp(-AUC(x=c(0, 1), y=c(0.207, 0.207), method="trapezoid"))))
  }
}

threshold.lambdas <- data.frame(VL = VLs, 
                             lambda = c(H1.threshold.lambda, H3.threshold.lambda), 
                             Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

threshold.lambda.plot <- ggplot(threshold.lambdas, aes(x=VL, y=lambda, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 10)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Force of infection") +
  theme_classic()

threshold.probs <- data.frame(VL = VLs, 
                           prob = c(H1.threshold.probs, H3.threshold.probs), 
                           lower = c(H1.threshold.lower, H3.threshold.lower), 
                           upper = c(H1.threshold.upper, H3.threshold.upper),
                           Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

threshold.prob.plot <- ggplot(threshold.probs, aes(x=VL, y=prob, color=Virus)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Virus, alpha=0.3)) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  guides(color="none", fill="none", alpha="none") +
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()


# HILL --------------------------------------------------------------------

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## range of parameter values to try
q_vals <- seq(1, 10, 0.1)
ka_vals <- seq(1, 10, 0.1)
n_vals <- seq(1, 10, 0.1)

## log prob tracker
H1N1.log.probs <- array(data=NA, dim=c(length(q_vals), length(ka_vals), length(n_vals)))

## nested loop
for (q in q_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H1N1_donor_names)))
  names(tmp.probs) <- H1N1_donor_names
  
  for (ka in ka_vals){
    for (n in n_vals){
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
            tmp.below_LOD <- which(tmp.yvals < LOD)
            if (length(tmp.below_LOD) >= 1){
              tmp.yvals[tmp.below_LOD] <- 0
            }
            ## Hill function
            tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
            tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.theta*q, method="trapezoid")
          }
          ## find integral from last neg to first positive test
          tmp.df <- H1N1.ferret.preds[[ferret]] %>%
            filter(time >= tmp.neg.time) %>%
            filter(time <= tmp.pos.time)
          tmp.yvals <- tmp.df$titer
          ## when titer is at or below LOD, AUC is 0
          tmp.below_LOD <- which(tmp.yvals < LOD)
          if (length(tmp.below_LOD) >= 1){
            tmp.yvals[tmp.below_LOD] <- 0
          }
          tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
          tmp.AUC <- AUC(tmp.df$time, tmp.theta*q, method="trapezoid")
          ## find prob of getting infected between last neg and first pos test
          tmp.interval.prob <- 1 - exp(-tmp.AUC)
          ## overall prob is 
          ## prob of not getting infected between 1dpi and last negative test 
          ## TIMES interval prob
          tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
        } else { ## if contact is uninfected, use this to calc prob of s
          tmp.yvals <- H1N1.ferret.preds[[ferret]]$titer
          ## when titer is at or below LOD, AUC is 0
          tmp.yvals[which(H1N1.ferret.preds[[ferret]]$titer < LOD)] <- 0
          tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
          tmp.AUC <- AUC(H1N1.ferret.preds[[ferret]]$time, tmp.theta*q, method="trapezoid")
          tmp.prob <- exp(-tmp.AUC)
        }
        tmp.probs[ferret] <- tmp.prob
      }
      ## replace any 0 prob with a very small number
      tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
      tmp.sum.prs <- sum(log(tmp.probs))
      ## for each q value, record likelihood of all ka and n values
      H1N1.log.probs[[which(q_vals == q), which(ka_vals == ka), which(n_vals == n)]] <- tmp.sum.prs
    }
    }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", q, "done"))
}  

## log prob tracker
H3N2.log.probs <- array(data=NA, dim=c(length(q_vals), length(ka_vals), length(n_vals)))

## nested loop
for (q in q_vals){
  ## create vector to track prob of transmission
  tmp.probs <- c(rep(NA, length(H3N2_donor_names)))
  names(tmp.probs) <- H3N2_donor_names
  
  for (ka in ka_vals){
    for (n in n_vals){
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
            tmp.below_LOD <- which(tmp.yvals < LOD)
            if (length(tmp.below_LOD) >= 1){
              tmp.yvals[tmp.below_LOD] <- 0
            }
            ## Hill function
            tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
            tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.theta*q, method="trapezoid")
          }
          ## find integral from last neg to first positive test
          tmp.df <- H3N2.ferret.preds[[ferret]] %>%
            filter(time >= tmp.neg.time) %>%
            filter(time <= tmp.pos.time)
          tmp.yvals <- tmp.df$titer
          ## when titer is at or below LOD, AUC is 0
          tmp.below_LOD <- which(tmp.yvals < LOD)
          if (length(tmp.below_LOD) >= 1){
            tmp.yvals[tmp.below_LOD] <- 0
          }
          tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
          tmp.AUC <- AUC(tmp.df$time, tmp.theta*q, method="trapezoid")
          ## find prob of getting infected between last neg and first pos test
          tmp.interval.prob <- 1 - exp(-tmp.AUC)
          ## overall prob is 
          ## prob of not getting infected between 1dpi and last negative test 
          ## TIMES interval prob
          tmp.prob <- exp(-tmp.neg.integral) * tmp.interval.prob
        } else { ## if contact is uninfected, use this to calc prob of s
          tmp.yvals <- H3N2.ferret.preds[[ferret]]$titer
          ## when titer is at or below LOD, AUC is 0
          tmp.yvals[which(H3N2.ferret.preds[[ferret]]$titer < LOD)] <- 0
          tmp.theta <- (tmp.yvals)^n / (ka^n + (tmp.yvals)^n)
          tmp.AUC <- AUC(H3N2.ferret.preds[[ferret]]$time, tmp.theta*q, method="trapezoid")
          tmp.prob <- exp(-tmp.AUC)
        }
        tmp.probs[ferret] <- tmp.prob
      }
      ## replace any 0 prob with a very small number
      tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
      tmp.sum.prs <- sum(log(tmp.probs))
      ## for each q value, record likelihood of all ka and n values
      H3N2.log.probs[[which(q_vals == q), which(ka_vals == ka), which(n_vals == n)]] <- tmp.sum.prs
    }
  }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", q, "done"))
}

## note that these index values are different than the code listed above
## because they were rerun for more precise estimates
which(H1N1.log.probs==max(H1N1.log.probs,na.rm=T), arr.ind=T)
H1N1.log.probs[7,17,111]
H1N1.MLE <- c("q" = q_vals[7], "ka" = ka_vals[17], "n" = n_vals[111])

which(H3N2.log.probs==max(H3N2.log.probs,na.rm=T), arr.ind=T)
H3N2.log.probs[3,51,2]
H3N2.MLE <- c("q" = q_vals[3], "ka" = ka_vals[51], "n"= n_vals[2])

## plots

H1.hill.lambda <- c()
H1.hill.probs <- c()
H1.hill.lower <- c()
H1.hill.upper <- c()
H3.hill.lambda <- c()
H3.hill.lower <- c()
H3.hill.upper <- c()
H3.hill.probs <- c()

for (v in VLs){
  if (v < LOD){
    H1.hill.lambda <- append(H1.hill.lambda, 0)
    H1.hill.probs <- append(H1.hill.probs, 0)
    H1.hill.lower <- append(H1.hill.lower, 0)
    H1.hill.upper <- append(H1.hill.upper, 0)
    H3.hill.lambda <- append(H3.hill.lambda, 0)
    H3.hill.probs <- append(H3.hill.probs, 0)
    H3.hill.lower <- append(H3.hill.lower, 0)
    H3.hill.upper <- append(H3.hill.upper, 0)
  } else {
    H1.theta <- v^17 / ((3.6^17) + v^17)
    H1.lower.theta <- v^6 / ((2.3)^6 + v^6)
    H1.upper.theta <- v^25 / ((4.2)^25 +v^25)
    H1.hill.lambda <- append(H1.hill.lambda, AUC(x=c(0, 1), y=c(0.7*H1.theta, 0.7*H1.theta), method="trapezoid"))
    H1.hill.probs <- append(H1.hill.probs, 1 - exp(-AUC(x=c(0, 1), y=c(0.7*H1.theta, 0.7*H1.theta), method="trapezoid")))
    H1.hill.lower <- append(H1.hill.lower, 1 - exp(-AUC(x=c(0, 1), y=c(0.3*H1.lower.theta, 0.3*H1.lower.theta), method="trapezoid")))
    H1.hill.upper <- append(H1.hill.upper, 1 - exp(-AUC(x=c(0, 1), y=c(1*H1.upper.theta, 1*H1.upper.theta), method="trapezoid")))
    
    H3.theta <- v^0.2 / ((13^0.2) + v^0.2)
    H3.lower.theta <- v^0.01 / ((3^0.01) + v^0.01)
    H3.upper.theta <- v^0.7 / ((23^0.7) + v^0.7)
    H3.hill.lambda <- append(H3.hill.lambda, AUC(x=c(0, 1), y=c(0.3*H3.theta, 0.3*H3.theta), method="trapezoid"))
    H3.hill.probs <- append(H3.hill.probs, 1 - exp(-AUC(x=c(0, 1), y=c(0.3*H3.theta, 0.3*H3.theta), method="trapezoid")))
    H3.hill.lower <- append(H3.hill.lower, 1 - exp(-AUC(x=c(0, 1), y=c(0.2*H3.lower.theta, 0.2*H3.lower.theta), method="trapezoid")))
    H3.hill.upper <- append(H3.hill.upper, 1 - exp(-AUC(x=c(0, 1), y=c(0.5*H3.upper.theta, 0.5*H3.upper.theta), method="trapezoid")))
  }
}

H1.retwisted.lower <- c(H1.hill.lower[which(H1.hill.lower == H1.hill.upper)],
                        H1.hill.upper[which(H1.hill.lower > H1.hill.upper)],
                        H1.hill.lower[which(H1.hill.lower < H1.hill.upper)])
H1.retwisted.upper <- c(H1.hill.lower[which(H1.hill.lower == H1.hill.upper)],
                       H1.hill.lower[which(H1.hill.lower > H1.hill.upper)],
                       H1.hill.upper[which(H1.hill.lower < H1.hill.upper)])

hill.lambdas <- data.frame(VL = VLs, 
                                lambda = c(H1.hill.lambda, H3.hill.lambda), 
                                Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))


hill.lambda.plot <- ggplot(hill.lambdas, aes(x=VL, y=lambda, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 10)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Force of infection") +
  theme_classic()

hill.probs <- data.frame(VL = VLs, 
                         prob = c(H1.hill.probs, H3.hill.probs),
                         lower = c(H1.hill.lower, H3.hill.lower),
                         upper = c(H1.hill.upper, H3.hill.upper),
                          Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

hill.prob.plot <- ggplot(hill.probs, aes(x=VL, y=prob, color=Virus)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Virus, alpha=0.3)) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  guides(color="none", fill="none", alpha="none") +
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()


# all plots ---------------------------------------------------------------

H1.log.lambda <- c()
H1.log.probs <- c()
H1.log.lower <- c()
H1.log.upper <- c()
H3.log.lambda <- c()
H3.log.probs <- c()
H3.log.lower <- c()
H3.log.upper <- c()

for (v in VLs){
  if (v < LOD){
    H1.log.lambda <- append(H1.log.lambda, 0)
    H3.log.lambda <- append(H3.log.lambda, 0)
    H1.log.probs <- append(H1.log.probs, 0)
    H1.log.lower <- append(H1.log.lower, 0)
    H1.log.upper <- append(H1.log.upper, 0)
    H3.log.probs <- append(H3.log.probs, 0)
    H3.log.lower <- append(H3.log.lower, 0)
    H3.log.upper <- append(H3.log.upper, 0)
  } else {
    H1.log.lambda <- append(H1.log.lambda, AUC(x=c(0,1), y=c(v*0.111, v*0.111), method="trapezoid"))
    H3.log.lambda <- append(H3.log.lambda, AUC(x=c(0,1), y=c(v*0.047, v*0.047), method="trapezoid"))
    
    H1.log.probs <- append(H1.log.probs, (1 - exp(-AUC(x=c(0,1), y=c(v*0.111, v*0.111), method="trapezoid"))))
    H1.log.lower <- append(H1.log.lower, (1 - exp(-AUC(x=c(0,1), y=c(v*0.068, v*0.068), method="trapezoid"))))
    H1.log.upper <- append(H1.log.upper, (1 - exp(-AUC(x=c(0,1), y=c(v*0.172, v*0.172), method="trapezoid"))))
    H3.log.probs <- append(H3.log.probs, (1 - exp(-AUC(x=c(0,1), y=c(v*0.047, v*0.047), method="trapezoid"))))
    H3.log.lower <- append(H3.log.lower, (1 - exp(-AUC(x=c(0,1), y=c(v*0.024, v*0.024), method="trapezoid"))))
    H3.log.upper <- append(H3.log.upper, (1 - exp(-AUC(x=c(0,1), y=c(v*0.082, v*0.082), method="trapezoid"))))
  }
}

log.lambdas <- data.frame(VL = VLs, 
                             lambda = c(H1.log.lambda, H3.log.lambda), 
                             Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

log.lambda.plot <- ggplot(log.lambdas, aes(x=VL, y=lambda, color=Virus)) +
  geom_line(linewidth=2) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 10)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Force of infection") +
  theme_classic()

log.probs <- data.frame(VL = VLs, 
                           prob = c(H1.log.probs, H3.log.probs), 
                        lower = c(H1.log.lower, H3.log.lower), 
                        upper = c(H1.log.upper, H3.log.upper),
                           Virus = c(rep("H1N1", length(VLs)), rep("H3N2", length(VLs))))

log.prob.plot <- ggplot(log.probs, aes(x=VL, y=prob, color=Virus)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Virus, alpha=0.3)) +
  scale_color_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  scale_fill_manual(values = c(plot_colors[[1]], plot_colors[[2]]), labels=c("Cal/2009", "Hong Kong/1968")) + 
  guides(color="none", fill="none", alpha="none") +
  scale_x_continuous(limits=c(0, 8), breaks = seq(0, 8, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8))) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()

top <- ggarrange(log.lambda.plot, linear.lambda.plot, threshold.lambda.plot, hill.lambda.plot, ncol=4, 
          common.legend = T, legend = "bottom", labels=c("A", "B", "C", "D"))
bottom <- ggarrange(log.prob.plot, linear.prob.plot, threshold.prob.plot, hill.prob.plot, ncol=4, 
          labels = c("E", "F", "G", "H"))

ggarrange(top, bottom, nrow=2, common.legend = T, legend="bottom")
