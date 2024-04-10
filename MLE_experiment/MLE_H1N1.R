library(tidyverse)
library(DescTools)

## import data
ferrets <- read_csv("/home/esomsen/within-host/H1N1_raw_titer_data.csv", col_names = T, show_col_types = F)
colnames(ferrets) <- c("Ferret_ID", "DI_RC", "DI_RC_Pair", "Dose", "time", "nw_titer", "donor_dose")

DI_ferrets <- ferrets %>%
  ## keep only donor ferrets
  filter(DI_RC == "DI") %>%
  dplyr::select(Ferret_ID, Dose, time, nw_titer) %>%
  mutate(Ferret_ID = as.factor(Ferret_ID)) %>%
  ## create column for days post exposure
  mutate(dpe = time - 1)
donor_names <- unique(DI_ferrets$Ferret_ID)

LOD <- 0.5

## function to add linear interpolations between measured datapoints
interpolation <- function(row1, row2, data, interval){
  index_1 <- data[row1, ]
  index_2 <- data[row2, ]
  times <- seq(index_1$dpe, index_2$dpe, interval)
  preds <- seq(index_1$nw_titer, index_2$nw_titer, length.out=length(times))
  df <- data.frame(dpe = times,
                   nw_titer = preds)
  return (df)
}

interpolation_interval <- 0.1

ferret_preds <- vector("list", length(donor_names))
names(ferret_preds) <- donor_names

## generate predictions for each ferret
for (ferret in donor_names){
  ferret_data <- DI_ferrets %>%
    filter(Ferret_ID == ferret) %>%
    select(c(dpe, nw_titer))
  df_1.3 <- interpolation(1, 2, ferret_data, interpolation_interval)
  df_3.5 <- interpolation(2, 3, ferret_data, interpolation_interval)
  df_5.7 <- interpolation(3, 4, ferret_data, interpolation_interval)
  df_7.9 <- interpolation(4, 5, ferret_data, interpolation_interval)
  df_9.11 <- interpolation(5, 6, ferret_data, interpolation_interval)
  combo <- rbind(ferret_data, df_1.3, df_3.5, df_5.7, df_7.9, df_9.11)
  combo <- combo %>%
    ## get rid of duplicate rows
    distinct() %>%
    ## specify if values are interpolated or measured
    mutate(type = if_else(dpe %in% c(0, 2, 4, 6, 8, 10), "measured", "predicted")) %>%
    arrange(dpe) %>%
    ## ensure that times are numeric for future integration
    mutate(dpe = as.numeric(dpe))
  ferret_preds[[ferret]] <- combo
}

## create function for F(t)
calculate_pr_contact_pos <- function(lambda_integral){
  prob <- 1 - exp(-lambda_integral)
  return (prob)
}

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## range of s values to try
s_vals <- seq(0, 0.2, 0.001)

## log prob tracker
log_probs <- c()

for (i in s_vals){
  ##loop through each ferret
  for (ferret in donor_names){
    ferret_data <- ferret_preds[[ferret]]
    ## empty vector for integral values
    integral_values <- c()
    ## loop through each timepoint
    for (t in times){
      VL <- ferret_data %>%
        filter(dpe <= t)
      xval <- VL$dpe
      yval <- VL$nw_titer * i
      ## find which elements have titer <= 0.5
      below_LOD <- which(VL$nw_titer <= 0.5)
      ## if any titers below, set the values of those to 0
      if (length(below_LOD) >= 1){
        for (k in below_LOD){
          yval[k] <- 0
        }
      }
      foi_val <- AUC(xval, yval, method = "trapezoid")
      integral_values <- append(integral_values, foi_val)
    }
    ferret_data <- ferret_data %>%
      mutate(foi_log = integral_values) %>%
      mutate(pr_contact_pos_log = calculate_pr_contact_pos(foi_log))
    ferret_preds[[ferret]] <- ferret_data
  }
  
  ## create vector
  interval_prs <- c(rep(0, length(donor_names)))
  names(interval_prs) <- donor_names
  
  for (ferret in donor_names){
    df <- ferret_preds[[ferret]]
    partner <- ferrets %>%
      filter(DI_RC == "RC") %>%
      filter(DI_RC_Pair == ferret) %>%
      mutate(dpe = time - 1)
    if (max(partner$nw_titer) > 0.5){
      first_positive_index <- which.max(partner$nw_titer > LOD)
      first_positive_time <- as.numeric(as.character(partner[first_positive_index, "dpe"]))
      ## if the first positive test is at 2dpi, call the last negative time at 0
      ## by definition, the integral at time 0 is 0
      if (first_positive_index > 1){
        last_negative_index <- first_positive_index - 1
        last_negative_time <- as.numeric(as.character(partner[last_negative_index, "dpe"]))
        integral_vals <- df %>%
          filter(dpe %in% c(last_negative_time, first_positive_time))
        ## interval pr = Pr(not seroconverting by last recipient negative) * Pr(seroconverting between last recipient negative and first recipient positive)
        ## to calculate the integral from t1 to t2, use raw nasal wash titers times proposed s value
        yval <- integral_vals$nw_titer * i
        new_integral <- AUC(x = integral_vals$dpe, y=yval, method = "trapezoid")
        t1_t2_prob <- calculate_pr_contact_pos(new_integral)
        prob <- (1 - integral_vals[[1, "pr_contact_pos_log"]]) * t1_t2_prob
      } else {
        # last negative time is 0
        ## don't need to change integral here because t1>t2 is still 0>t2
        upper_integral_val <- df %>%
          filter(near(dpe, first_positive_time))
        prob <- upper_integral_val$pr_contact_pos_log * 1 ## 1-0, 
        ## because prob of not seroconverting by last recipient negative
        ## is 1 when the integral value is 0
      }
    } else { ## if contact is uninfected, use prob of contact remaining uninfected by end of study 
      prob <- 1 - df[length(df$dpe)-1, "pr_contact_pos_log"]
    }
    interval_prs[ferret] <- prob
  }
  
  interval_prs <- unlist(interval_prs)
  log_interval_prs <- log(interval_prs)
  ## replace probabilities of 0 with log of smallest number in R
  log_interval_prs <- replace(log_interval_prs, log_interval_prs==-Inf, -744.4401)
  interval_log_p <- sum(log_interval_prs)
  log_probs <- append(log_probs, interval_log_p)
  print(paste("Round", i, "done"))
}

H1N1_MLE_trace <- data.frame(s = s_vals, 
                           prob = log_probs)

max <- max(H1N1_MLE_trace$prob)
index <- match(max, H1N1_MLE_trace$prob)
MLE <- H1N1_MLE_trace[index,]

CI_cutoff <- max - 1.92
find_CIs <- near(CI_cutoff, H1N1_MLE_trace$prob, tol = 0.053) ## toggle tol to find the two values closest to the CI cutoff
CIs <- H1N1_MLE_trace[find_CIs,]

save.image(file = "s_MLE_H1N1.RData")
save(H1N1_MLE_trace, file="H1N1_MLE_trace.RData") 

## the below code generates the transmission probability over time for set values of s
## these results are displayed in fig 3b
## the required input is interpolated viral load data for each ferret
## which is generated above

## vector of times at which to calculate the integral
times <- seq(0, 10, 0.1)
## pick random value for s
s_vals <- c(0.05, 0.1, 0.2)

panel_3b <- vector("list", length(s_vals))

for (i in 1:length(panel_3b)){
  s <- s_vals[i]
  for (ferret in donor_names){
    ferret_data <- ferret_preds[[ferret]]
    ## empty vector for integral values
    integral_values <- c()
    ## loop through each timepoint
    for (t in times){
      VL <- ferret_data %>%
        filter(dpe <= t)
      xval <- VL$dpe
      yval <- VL$nw_titer * s
      ## find which elements have titer <= 0.5
      below_LOD <- which(VL$nw_titer <= 0.5)
      ## if any titers below, set the values of those to 0
      if (length(below_LOD) >= 1){
        for (k in below_LOD){
          yval[k] <- 0
        }
      }
      foi_val <- AUC(xval, yval, method = "trapezoid")
      integral_values <- append(integral_values, foi_val)
    }
    ferret_data <- ferret_data %>%
      mutate(foi_log = integral_values) %>%
      mutate(pr_contact_pos_log = calculate_pr_contact_pos(foi_log))
    ferret_preds[[ferret]] <- ferret_data
  }
  panel_3b[[i]] <- ferret_preds
}

save(panel_3b, file="panel_3b_data.RData") 
