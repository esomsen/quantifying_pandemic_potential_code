library(tidyverse)
library(DescTools)
library(ggpubr)

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
H3N2.MLE <- prob.trace[which.max(H3N2.log.probs),1]


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
H1N1.joint.log.probs[284,708]

which(H3N2.joint.log.probs == max(H3N2.joint.log.probs), arr.ind=T)
H3N2.MLE <- c("s" = s_vals[118], "h" = h_vals[1])
H3N2.joint.log.probs[1,118]

# HILL --------------------------------------------------------------------

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## range of parameter values to try
k_vals <- seq(1, 10, 0.1)
ka_vals <- seq(1, 10, 0.1)
n_vals <- seq(1, 10, 0.1)

## log prob tracker
H1N1.log.probs <- array(data=NA, dim=c(length(k_vals), length(ka_vals), length(n_vals)))

## nested loop
for (k in k_vals){
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
            tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.theta*k, method="trapezoid")
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
          tmp.AUC <- AUC(tmp.df$time, tmp.theta*k, method="trapezoid")
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
          tmp.AUC <- AUC(H1N1.ferret.preds[[ferret]]$time, tmp.theta*k, method="trapezoid")
          tmp.prob <- exp(-tmp.AUC)
        }
        tmp.probs[ferret] <- tmp.prob
      }
      ## replace any 0 prob with a very small number
      tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
      tmp.sum.prs <- sum(log(tmp.probs))
      ## for each k value, record likelihood of all ka and n values
      H1N1.log.probs[[which(k_vals == k), which(ka_vals == ka), which(n_vals == n)]] <- tmp.sum.prs
    }
    }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", k, "done"))
}  

## log prob tracker
H3N2.log.probs <- array(data=NA, dim=c(length(k_vals), length(ka_vals), length(n_vals)))

## nested loop
for (k in k_vals){
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
            tmp.neg.integral <- AUC(tmp.neg.titers$time, tmp.theta*k, method="trapezoid")
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
          tmp.AUC <- AUC(tmp.df$time, tmp.theta*k, method="trapezoid")
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
          tmp.AUC <- AUC(H3N2.ferret.preds[[ferret]]$time, tmp.theta*k, method="trapezoid")
          tmp.prob <- exp(-tmp.AUC)
        }
        tmp.probs[ferret] <- tmp.prob
      }
      ## replace any 0 prob with a very small number
      tmp.probs <- replace(tmp.probs, tmp.probs==0, 2.225074e-308)
      tmp.sum.prs <- sum(log(tmp.probs))
      ## for each k value, record likelihood of all ka and n values
      H3N2.log.probs[[which(k_vals == k), which(ka_vals == ka), which(n_vals == n)]] <- tmp.sum.prs
    }
  }
  rm(list=ls(pattern="^tmp"))
  print(paste("Round", k, "done"))
}


which(x==max(x,na.rm=T), arr.ind=T)

H1N1.MLE <- prob.trace[which.max(H1N1.log.probs),1]
H3N2.MLE <- prob.trace[which.max(H3N2.log.probs),1]



# EXAMPLE PLOTS -----------------------------------------------------------


VLs <- seq(0, 10, 0.1)
LOD <- 1

s.log <- 0.1
s.linear <- 0.000001

log.probs <- c()
log.lambda <- c()
linear.probs <- c()
linear.lambda <- c()

for (v in VLs){
  if (v < LOD){
    linear.probs <- append(linear.probs, 0)
    linear.lambda <- append(linear.lambda, 0)
    log.probs <- append(log.probs, 0)
    log.lambda <- append(log.lambda, 0)
  } else {
    linear.probs <- append(linear.probs, (1 - exp(-AUC(x=c(0,1), y=c((10^v)*s.linear, (10^v)*s.linear), method="trapezoid"))))
    linear.lambda <- append(linear.lambda, AUC(x=c(0,1), y=c((10^v)*s.linear, (10^v)*s.linear), method="trapezoid"))
    log.probs <- append(log.probs, (1 - exp(-AUC(x=c(0,1), y=c(v*s.log, v*s.log), method="trapezoid"))))
    log.lambda <- append(log.lambda, AUC(x=c(0,1), y=c(v*s.log, v*s.log), method="trapezoid"))
  }
}

threshold.probs <- c()
threshold.lambda <- c()
s.threshold <- 0.3
h <- 4

for (v in VLs){
  if (v < h){
    threshold.probs <- append(threshold.probs, 0)
    threshold.lambda <- append(threshold.lambda, 0)
  } else {
    threshold.probs <- append(threshold.probs, (1 - exp(-AUC(x=c(0, 1), y=c(s.threshold, s.threshold), method="trapezoid"))))
    threshold.lambda <- append(threshold.lambda, AUC(x=c(0, 1), y=c(s.threshold, s.threshold), method="trapezoid"))
  }
}

hill.probs <- c()
hill.lambda <- c()
k <- 5
ka <- 5
n <- 3

for (v in VLs){
  if (v < LOD){
    hill.probs <- append(hill.probs, 0)
    hill.lambda <- append(hill.lambda, 0)
  } else {
    theta <- v^n / ((ka^n) + v^n)
    hill.probs <- append(hill.probs, 1 - exp(-AUC(x=c(0, 1), y=c(k*theta, k*theta), method="trapezoid")))
    hill.lambda <- append(hill.lambda, AUC(x=c(0, 1), y=c(k*theta, k*theta), method="trapezoid"))
  }
}

linear.plot <- ggplot(data.frame(VL = VLs, 
                                 prob = linear.probs, 
                                 LOD = ifelse(VLs < LOD, "y", "n")), aes(x=VL, y=prob, color=LOD)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=c("black", "grey")) +
  guides(color="none") +
  scale_x_continuous(limits=c(0, 10), breaks = seq(0, 10, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8), expression(10^10))) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()

log.plot <- ggplot(data.frame(VL = VLs, 
                              prob = log.probs, 
                              LOD = ifelse(VLs < LOD, "y", "n")), aes(x=VL, y=prob, color=LOD)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=c("black", "grey")) +
  guides(color="none") +
  scale_x_continuous(limits=c(0, 10), breaks = seq(0, 10, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8), expression(10^10))) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()

threshold.plot <- ggplot(data.frame(VL = VLs, 
                                    prob = threshold.probs,
                                    LOD = ifelse(VLs < LOD, "y", "n")), aes(x=VL, y=prob, color=LOD)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=c("black", "grey")) +
  guides(color="none") +
  scale_x_continuous(limits=c(0, 10), breaks = seq(0, 10, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8), expression(10^10))) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()

hill.plot <- ggplot(data.frame(VL = VLs, 
                               prob = hill.probs,
                               LOD = ifelse(VLs < LOD, "y", "n")), aes(x=VL, y=prob, color=LOD)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=c("black", "grey")) +
  guides(color="none") +
  scale_x_continuous(limits=c(0, 10), breaks = seq(0, 10, 2), labels=c(expression(10^0), expression(10^2), expression(10^4), expression(10^6), expression(10^8), expression(10^10))) +
  labs(x=expression(paste("Viral titer (", TCID[50], "/mL)")), y="Probability of transmission") +
  theme_classic()

ggarrange(linear.plot, log.plot, threshold.plot, hill.plot, ncol=4, labels=c("A", "B", "C", "D"))

