
# packages
# library(boot)
# library(mice)
# library(dplyr)
# library(zoo)
# library(MASS)
# library(ggplot2)
# library(survival)

# Code for dynamic treatment regimes comparing P/F>150 vs >200 vs >250 depending on whether PEEP is <8 cmH2O or >=8cmH2O
# For primary outcome: Restricted mean time lost (RMTL) of successful extubation


num_iterations = 300
df_survival


count <- 0

fboot <- function(df, indices) {
  
  total_ids <- unique(df$stay_id) # get all distinct patients
  total_ids_df <- data.frame(stay_id = total_ids)
  sample_ids <- total_ids_df[indices,] # allows boot to select sample
  
  df <- df[df$stay_id %in% sample_ids,] 
  
  print(count)
  assign("count", count+1, envir = .GlobalEnv)
  
  
  
  df_after_TTE_period = df
  ################################################################################
  # clean up data for joining after weights are found in TTE period
  ###############################################################################
  
  # convert 'sedation' to factors
  df_after_TTE_period$sedation <- factor(df_after_TTE_period$sedation, levels = c('Awake', 'Moderate sedation', 'Deep sedation'), ordered = TRUE)
  df_after_TTE_period$sedation_w_missing <- factor(df_after_TTE_period$sedation_w_missing, levels = c('Awake', 'Moderate sedation', 'Deep sedation'), ordered = TRUE)
  
  # Make sure cv_sofa score in 'extubated' patients is latest known cv_sofa score
  df_after_TTE_period <- df_after_TTE_period %>%
    arrange(stay_id, hour) %>%  
    group_by(stay_id) %>%
    mutate(cv_sofa = na.locf(cv_sofa, na.rm = FALSE)) # latest known cv_sofa
  
  ################################################################################
  
  # convert 'cv_sofa' to factors
  df_after_TTE_period$cv_sofa <- cut(
    df_after_TTE_period$cv_sofa,
    breaks = c(-Inf, 1, 3, Inf),   # intervals: less than 1, 1 to 3, greater than 3
    labels = c("<1", "1-3", ">3"),
    right = FALSE,                 # intervals are left-inclusive, right-exclusive: [a,b)
    ordered_result = TRUE          # creates an ordered factor
  )
  
  df_after_TTE_period$cv_sofa_w_missing <- cut(
    df_after_TTE_period$cv_sofa_w_missing,
    breaks = c(-Inf, 1, 3, Inf),   # intervals: less than 1, 1 to 3, greater than 3
    labels = c("<1", "1-3", ">3"),
    right = FALSE,                 # intervals are left-inclusive, right-exclusive: [a,b)
    ordered_result = TRUE          # creates an ordered factor
  )
  
  ################################################################################
  # Find TTE period: Q3 + 1.5 * IQR
  ################################################################################
  # find hour of first switch (first find eligible time)
  df <- df %>%
    group_by(stay_id) %>%
    mutate(
      eligible_time = if (any(hour == 0)) time[which(hour == 0)[1]] else as.POSIXct(NA)
    ) %>%
    ungroup()
  
  df <- df %>%
    mutate(
      hour_of_switch = as.numeric(first_switch_time - eligible_time, units = "hours")
    )
  
  
  df_summary <- df %>%
    group_by(stay_id) %>%
    summarise(hour_of_switch = first(hour_of_switch), .groups = "drop")
  
  Q1 <- quantile(df_summary$hour_of_switch, 0.25, na.rm = TRUE)
  Q3 <- quantile(df_summary$hour_of_switch, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  upper_bound <- Q3 + 1.5 * IQR
  
  ################################################################################
  # Set TTE period 
  
  df <- df[-which(df$hour > upper_bound),]
  
  ###############################################################################
  # DEAL WITH MISSING DATA
  ###############################################################################
  
  # convert 'sedation' to factors
  df$sedation <- factor(df$sedation, levels = c('Awake', 'Moderate sedation', 'Deep sedation'), ordered = TRUE)
  df$sedation_w_missing <- factor(df$sedation_w_missing, levels = c('Awake', 'Moderate sedation', 'Deep sedation'), ordered = TRUE)
  
  # use missing imputation for all controlled/assisted rows only
  imv_data <- df[!((df$extubated == 1) | (df$died == 1)),]
  no_imv_data <- df[(df$extubated == 1) | (df$died == 1), ]
  
  # Use MICE for missingness imputation
  cols_to_impute <- c("age","PF_ratio", "ph","pco2", "peep", "driving", "rr_set","vt", "cv_sofa", "sedation", "nmb")
  imputed_data <- tryCatch({
    complete(mice(imv_data[, cols_to_impute], printFlag = FALSE))
  }, error = function(e) {
    warning("Error occurred during mice imputation. Switching to median imputation.")
    imputed_data_median <- imv_data
    for (col in cols_to_impute) {
      if (anyNA(imputed_data_median[[col]])) {
        if (is.factor(imputed_data_median[[col]])) {
          # Convert factor to numeric (assuming the levels can be coerced to numeric)
          imputed_data_median[[col]] <- as.numeric(as.character(imputed_data_median[[col]]))
        }
        imputed_data_median[[col]] <- ifelse(is.na(imputed_data_median[[col]]), median(imputed_data_median[[col]], na.rm = TRUE), imputed_data_median[[col]])
      }
    }
    # convert 'sedation' back to a factors
    imputed_data_median$sedation <- factor(imputed_data_median$sedation, ordered = TRUE)
    imputed_data_median$sedation_w_missing <- factor(imputed_data_median$sedation_w_missing, ordered = TRUE)
    
    return(imputed_data_median)
  })
  
  
  imv_data[, cols_to_impute] <- imputed_data
  
  # combine the imputed imv data and the no imv data
  df <- rbind(imv_data, no_imv_data)
  
  # Make sure cv_sofa score in 'extubated' patients is latest known cv_sofa score
  df <- df %>%
    arrange(stay_id, hour) %>%  
    group_by(stay_id) %>%
    mutate(cv_sofa = na.locf(cv_sofa, na.rm = FALSE)) # latest known cv_sofa
  
  ################################################################################
  
  # convert 'cv_sofa' to factors
  df$cv_sofa <- cut(
    df$cv_sofa,
    breaks = c(-Inf, 1, 3, Inf),   # intervals: less than 1, 1 to 3, greater than 3
    labels = c("<1", "1-3", ">3"),
    right = FALSE,                 # intervals are left-inclusive, right-exclusive: [a,b)
    ordered_result = TRUE          # creates an ordered factor
  )
  
  df$cv_sofa_w_missing <- cut(
    df$cv_sofa_w_missing,
    breaks = c(-Inf, 1, 3, Inf),   # intervals: less than 1, 1 to 3, greater than 3
    labels = c("<1", "1-3", ">3"),
    right = FALSE,                 # intervals are left-inclusive, right-exclusive: [a,b)
    ordered_result = TRUE          # creates an ordered factor
  )
  
  
  ################################################################################
  
  # add columns of PF, driving pressure, ph, peep values and still_controlled on previous time point (time-1)
  df <- df %>%
    arrange(stay_id, hour) %>%
    group_by(stay_id) %>%
    mutate(PF_ratio_timemin1 = lag(PF_ratio, default = NA),
           driving_timemin1  = lag(driving, default = NA),
           vt_timemin1  = lag(vt, default = NA),
           pco2_timemin1  = lag(pco2, default = NA),
           rr_set_timemin1  = lag(rr_set, default = NA),
           ph_timemin1  = lag(ph, default = NA), 
           peep_timemin1  = lag(peep, default = NA),
           cv_sofa_timemin1  = lag(cv_sofa, default = NA),
           sedation_timemin1  = lag(sedation, default = NA),
           switched_timemin1  = lag(switched, default = NA),
           nmb_timemin1  = lag(nmb, default = NA))
  
  df$switched_timemin1 <- ifelse(is.na(df$switched_timemin1), 0, df$switched_timemin1)
  
  # fill ph, peep, PF, dp , min_vol for extubated/death with latest known value, for confounding
  df <- df %>%
    arrange(stay_id, hour) %>%
    group_by(stay_id) %>%
    mutate(ph_timemin1 = zoo::na.locf(ph_timemin1, na.rm = FALSE)) %>%
    mutate(pco2_timemin1 = zoo::na.locf(pco2_timemin1, na.rm = FALSE)) %>%
    mutate(peep_timemin1 = zoo::na.locf(peep_timemin1, na.rm = FALSE)) %>%
    mutate(PF_ratio_timemin1 = zoo::na.locf(PF_ratio_timemin1, na.rm = FALSE)) %>%
    mutate(driving_timemin1 = zoo::na.locf(driving_timemin1, na.rm = FALSE)) %>%
    mutate(rr_set_timemin1 = zoo::na.locf(rr_set_timemin1, na.rm = FALSE)) %>%
    mutate(vt_timemin1 = zoo::na.locf(vt_timemin1, na.rm = FALSE)) %>%
    mutate(sedation_timemin1 = zoo::na.locf(sedation_timemin1, na.rm = FALSE))%>%
    mutate(nmb_timemin1 = zoo::na.locf(nmb_timemin1, na.rm = FALSE))
  
  
  
  ################################################################################
  # Propensity score model
  ################################################################################
  print('model')
  
  
  data_model <- df[df$switched_timemin1 == 0 & df$hour != 0 ,]
  
  # find probability of switching
  ps_model <- glm(formula = switched ~ factor(hour) + ph_timemin1*pco2_timemin1 + pco2_timemin1*factor(hour) + ph_timemin1*factor(hour) + PF_ratio_timemin1*factor(hour)+ peep_timemin1*factor(hour) + driving_timemin1*factor(hour) + rr_set_timemin1*factor(hour) + vt_timemin1*factor(hour)+ cv_sofa_timemin1*factor(hour) + sedation_timemin1*factor(hour),
                  data = data_model, # exclude day 0 as there are no predictor variables (no hour before)
                  family = "binomial")
  
  
  df$ps <- NA # initialize column
  df$ps[df$switched_timemin1 == 0 & df$hour != 0 ] <- predict(ps_model,type = "response")
  df$ps[df$switched_timemin1 == 1] <- 1
  df$ps[df$hour == 0] <- 0 # probability of switching at day=0 is 0
  
  
  ################################################################################
  # Define considered regime *thresholds*
  ################################################################################
  
  
  
  threshold_values <- list(
    X_lowPEEP = c(150, 200, 250), # smaller than
    X_highPEEP = c(150, 200, 250) # larger than
  )
  
  # Generate all possible combinations of threshold values
  combinations <- expand.grid(
    X_lowPEEP = threshold_values$X_lowPEEP,
    X_highPEEP = threshold_values$X_highPEEP
    
  )
  # Add a count column 'regime' starting at 1 for the first combination
  combinations <- cbind(regime = rownames(combinations), combinations)
  
  # Create vectors storing information for each regime
  percentage_uncensored_values <- numeric(nrow(combinations))
  cuminc_suc_ext_values <- numeric(nrow(combinations))
  RMTL_suc_ext_values <- numeric(nrow(combinations))
  PF_lowPEEP_values <- numeric(nrow(combinations))
  PF_highPEEP_values <- numeric(nrow(combinations))
  switch_failed_values <- numeric(nrow(combinations))
  
  
  # store results of each time point for plotting
  ci_df <- as.data.frame(matrix(NA, nrow = length(threshold_values$X_lowPEEP) * length(threshold_values$X_highPEEP), ncol = length(seq(0, 672, by = 24))))
  
  
  #################################################################################
  
  
  combinations_sorted <- combinations[order(as.integer(combinations$regime)), ]
  
  # Loop over each combination
  for (i in rev(1:nrow(combinations_sorted))) {
  
    regime_nr <- combinations_sorted$regime[i]
    
    # Set the threshold values based on the current combination
    X_lowPEEP <- combinations_sorted$X_lowPEEP[i]
    X_highPEEP <- combinations_sorted$X_highPEEP[i]
    
    print(paste('regime:', regime_nr, " PEEP<8, PF>", X_lowPEEP, "; PEEP>8, PF>", X_highPEEP))
    
    # create a copy of the dataset and assign it to regime 1
    df_csw <- df
    # create new column specifying which regime
    df_csw$regime <- regime_nr
    
    ################################################################################
    # Censor patients the moment they are no longer compatible with regime
    # and delete rows after being censored
    ################################################################################
    
    
    # add columns eligible_timemin1
    # eligible if: PEEP<8 & P/F>XlowPEEP, or PEEP >=8 and P/F>XhighPEEP
    df_csw <- df_csw %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      # eligible on time-1?
      mutate(eligible_timemin1 = ifelse(
        (switched_timemin1 != 1 ) & (nmb_timemin1 != 1) & 
          (
            (!is.na(peep_timemin1) & (peep_timemin1 < 8) & !is.na(PF_ratio_timemin1) & (PF_ratio_timemin1 > X_lowPEEP)) |
              (!is.na(peep_timemin1) & (peep_timemin1 >= 8) & !is.na(PF_ratio_timemin1) & (PF_ratio_timemin1 > X_highPEEP))
          ),
        1, 0)
      )
    # censor if eligible but did not switch, or switched but not eligible
    df_csw <- df_csw %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      mutate(censored = ifelse(
        ((switched_timemin1 !=1) & (switched == 1) & (eligible_timemin1 == 0)) | # sc1 – switches but not eligible for switching 
          ((switched != 1) & (died != 1) & (ltfu_discharged != 1) & (eligible_timemin1 == 1)),  # sc4 – does not switch but eligible 
        1, 0)
      )  
    
    
    # remove all rows after the first censoring occurrence
    df_csw <- df_csw %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      mutate(before_censoring_occurrence = cumsum(censored == 1) ==0)%>% # before censoring occurs --> TRUE
      mutate(first_censoring_occurrence = cumsum(censored == 1) == 1) %>% # first censoring occurrence --> TRUE
      filter(before_censoring_occurrence | (first_censoring_occurrence & (censored == 1))) %>%
      dplyr::select(-c(before_censoring_occurrence, first_censoring_occurrence))
    
    # find probability of remaining uncensored for stabilization
    df_csw$uncensor <- 1-df_csw$censored
    cens_model <- glm(formula = uncensor ~ factor(hour) , 
                      data = df_csw,  
                      family = "binomial")
    df_csw$prob_uncen <- NA # initialize column
    df_csw$prob_uncen <- predict(cens_model, df_csw, type = "response")
    
    # compute prob of remaining uncensored
    df_csw <- df_csw %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      mutate(
        prob_uncen3 = prob_uncen^uncensor * (1-prob_uncen)^(1-uncensor)) # when censored, 1-prob uncensored
    
    ##############################################################################
    # find average PEEP and PF just before switch
    ##############################################################################
    df_switch <- df_csw[(df_csw$censored == 0) & (df_csw$switched == 1) & (df_csw$switched_timemin1 != 1), ]
    
    # Split into two datasets
    df_lowPEEP  <- df_switch[df_switch$peep_timemin1 <8 ,] 
    df_highPEEP <- df_switch[df_switch$peep_timemin1 >=8 ,]
    
    median_val <- median(df_lowPEEP$PF_ratio_timemin1, na.rm = TRUE)
    Q1 <- quantile(df_lowPEEP$PF_ratio_timemin1, probs = c(0.25), na.rm = TRUE)
    Q3 <- quantile(df_lowPEEP$PF_ratio_timemin1, probs = c(0.75), na.rm = TRUE)
    PFspread_lowPEEP <- paste0(round(median_val), " (", round(Q1), "–", round(Q3), ")")
    
    median_val <- median(df_highPEEP$PF_ratio_timemin1, na.rm = TRUE)
    Q1 <- quantile(df_highPEEP$PF_ratio_timemin1, probs = c(0.25), na.rm = TRUE)
    Q3 <- quantile(df_highPEEP$PF_ratio_timemin1, probs = c(0.75), na.rm = TRUE)
    PFspread_highPEEP <- paste0(round(median_val), " (", round(Q1), "–", round(Q3), ")")
    
    
    ################################################################################
    ### IPCW
    ################################################################################
    
    # compute prob of remaining uncensored
    df_csw <- df_csw %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      mutate(
        prob_uncen2 = ps^switched * (1-ps)^(1-switched))
    
    # Compute (stabilized) Inverse Probability of Compatibility Weights using propensity scores and decision rules for different regimes. 
    df_csw <- mutate(group_by(df_csw, stay_id), w_IPW = cumprod(prob_uncen3) / cumprod(prob_uncen2))
    
    # truncate weights such that all weights >10 are set to 10
    trunc.cutoff <- 10
    df_csw$w_IPW <- ifelse(df_csw$w_IPW>trunc.cutoff, trunc.cutoff, df_csw$w_IPW)
    
    
    ##############################################################################
    # Combine with data after the TTE period and make sure censored and weights are carried forward
    
    colnames(df_csw) 
    colnames(df_after_TTE_period)
    
    # Get the stay_ids already in df_csw
    stay_ids <- unique(df_csw$stay_id)
    
    # Extract rows from df_survival beyond TTE_period
    df_post_TTE_period <- df_after_TTE_period %>%
      filter(stay_id %in% stay_ids) %>%
      dplyr::anti_join(df_csw, by = c("stay_id", "hour"))
    
    
    # Combine the TTE period data with extended rows
    df_combined <- bind_rows(df_csw, df_post_TTE_period) %>%
      arrange(stay_id, hour)
    
    
    # Apply LOCF to 'censored' and 'weight' by stay_id
    df_combined <- df_combined %>%
      arrange(stay_id, hour) %>%
      group_by(stay_id) %>%
      tidyr::fill(censored, w_IPW, .direction = "down") %>%
      ungroup()
    
    summary(df_combined)
    
    
    ################################################################################
    ### Aalen-Johanson estimator
    ################################################################################
    
    
    percentage_uncensored <- df_combined %>%
      group_by(stay_id) %>%
      filter(all(censored == 0)) %>%
      summarise() %>%
      summarise(pct = n() / n_distinct(df_combined$stay_id) * 100) %>%
      pull(pct)
    
    # Keep only non-censored patients
    df_combined <- df_combined[df_combined$censored == 0,] 
    
    # among those who switched, find percentage of failed
    switch_failed_perc <- df_combined %>%
      group_by(stay_id) %>%
      summarise(
        switched = any(switched == 1),
        switch_failed = any(switch_failed == 1)
      ) %>%
      filter(switched) %>%
      summarise(pct = mean(switch_failed) * 100) %>%
      pull(pct)
    
    df_combined <- df_combined %>%
      arrange(stay_id, hour) %>%  
      group_by(stay_id) %>%
      mutate(
        T_start = coalesce(lag(hour), -1),   # previous time point
        T_stop = hour
      ) %>%
      ungroup()
    
    df_combined$event <- df_combined$successfully_extubated * 1 + df_combined$died * 2
    df_combined$event <- factor(df_combined$event, 0:2, c("atrisk", "extubated_success", "mortality")) 
    
    
    ci <- survfit(Surv(T_start, T_stop, event) ~ 1, data = df_combined, weights = w_IPW, id = stay_id)
    
    
    ci_28 <- summary(ci, times=28*24, extend=TRUE)$pstate[,2]
    # find RMTL
    rmtl_28 <- (summary(ci, rmean=28*24)$table[2,3])/24
    
    # Store results at each time point for plot
    time_points_plot <- seq(0, 672, by = 24)
    ci_28_plot <- summary(ci, times=time_points_plot, extend=TRUE)$pstate[,2]
    ci_df[i, 1:length(ci_28_plot)] <- ci_28_plot
    
    
    percentage_uncensored_values[i] <- percentage_uncensored
    cuminc_suc_ext_values[i] <- ci_28
    RMTL_suc_ext_values[i] <- rmtl_28
    PF_lowPEEP_values[i] <- PFspread_lowPEEP
    PF_highPEEP_values[i] <- PFspread_highPEEP
    switch_failed_values[i] <- switch_failed_perc
    
    
    
  }
  
  res <- c(percentage_uncensored_values, cuminc_suc_ext_values, RMTL_suc_ext_values, switch_failed_values, unlist(ci_df[ , paste0("V", 1:29)]))
  names(res) <- c("NR_stays_perc", "cuminc_suc_ext","RMTL_suc_ext", "perc_switch_fail", colnames(ci_df))
  return(res)
}

# call the function
results_boot_primaryoutcome <- boot(data=df_survival, statistic=fboot, R=num_iterations)
save.image(file = file.path(output_folder_name,"/workspace_DTR_primary.RData"))

# create combinations df
threshold_values <- list(
  X_lowPEEP = c(150, 200, 250), # larger than
  X_highPEEP = c(150, 200, 250) # larger than
)
combinations <- expand.grid(
  X_lowPEEP = threshold_values$X_lowPEEP,
  X_highPEEP = threshold_values$X_highPEEP)
combinations <- cbind(regime = rownames(combinations), combinations)

combinations_sorted <- combinations[order(as.integer(combinations$regime)), ]


# find colmeans
nr_regimes = nrow(combinations_sorted)
nr_stays_perc <- colMeans(results_boot_primaryoutcome$t)[1:nr_regimes]
cuminc_suc_ext <- colMeans(results_boot_primaryoutcome$t)[(nr_regimes+1):(nr_regimes*2)]
RMTL_suc_ext <- colMeans(results_boot_primaryoutcome$t)[(nr_regimes*2+1):(nr_regimes*3)]
perc_switch_fail <- colMeans(results_boot_primaryoutcome$t)[(nr_regimes*3+1):(nr_regimes*4)]

# find CIs:

# Initialize list to store formatted strings
ci_strings_nr_stays_perc <- vector("character", length = nr_regimes)
# Loop over each index and compute confidence intervals
for (i in 1:nr_regimes) {
  ci <- boot.ci(results_boot_primaryoutcome, type="perc", index=i)$percent[4:5]
  
  mean_val <- round(colMeans(results_boot_primaryoutcome$t)[i], 3)
  lower_ci <- round(ci[1], 1)
  upper_ci <- round(ci[2], 1)
  
  print(paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")"))
  
  ci_strings_nr_stays_perc[i] <- paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")")
}

# Initialize list to store formatted strings
ci_strings_cuminc_suc_ext <- vector("character", length = nr_regimes)
# Loop over each index and compute confidence intervals
for (i in (nr_regimes+1):(nr_regimes*2)) {
  ci <- boot.ci(results_boot_primaryoutcome, type="perc", index=i)$percent[4:5]
  
  mean_val <- round(colMeans(results_boot_primaryoutcome$t)[i] *100, 1)
  lower_ci <- round(ci[1] *100, 1)
  upper_ci <- round(ci[2] *100, 1)
  
  print(paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")"))
  
  ci_strings_cuminc_suc_ext[i-nr_regimes] <- paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")")
}

# Initialize list to store formatted strings
ci_strings_RMTL_suc_ext <- vector("character", length = nr_regimes)
# Loop over each index and compute confidence intervals
for (i in (nr_regimes*2+1):(nr_regimes*3)) {
  ci <- boot.ci(results_boot_primaryoutcome, type="perc", index=i)$percent[4:5]
  
  mean_val <- round(colMeans(results_boot_primaryoutcome$t)[i], 1)
  lower_ci <- round(ci[1], 1)
  upper_ci <- round(ci[2], 1)
  
  print(paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")"))
  
  ci_strings_RMTL_suc_ext[i-(nr_regimes*2)] <- paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")")
}

# Initialize list to store formatted strings
ci_strings_switch_fail_perc <- vector("character", length = nr_regimes)
# Loop over each index and compute confidence intervals
for (i in (nr_regimes*3+1):(nr_regimes*4)) {
  ci <- boot.ci(results_boot_primaryoutcome, type="perc", index=i)$percent[4:5]
  
  mean_val <- round(colMeans(results_boot_primaryoutcome$t)[i], 1)
  lower_ci <- round(ci[1], 1)
  upper_ci <- round(ci[2], 1)
  
  print(paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")"))
  
  ci_strings_switch_fail_perc[i-(nr_regimes*3)] <- paste0(mean_val, " (", lower_ci, " – ", upper_ci, ")")
}

combinations_sorted$boot_nr_stays_perc <- ci_strings_nr_stays_perc
combinations_sorted$boot_cuminc_suc_ext <- ci_strings_cuminc_suc_ext
combinations_sorted$boot_RMTL_suc_ext <- ci_strings_RMTL_suc_ext
combinations_sorted$boot_switch_fail_perc <- ci_strings_switch_fail_perc

write.xlsx(combinations_sorted, file.path(output_folder_name,"/boot_primary.xlsx"))


################################################################################
# plot cumulative incidences

# subset results for plotting
col_id_plot_start <- nr_regimes * 4 + 1
boot_results_plot <- results_boot_primaryoutcome$t[, col_id_plot_start:ncol(results_boot_primaryoutcome$t)]

nr_timepoints <- 29
nr_regimes    <- nr_regimes
nr_boot       <- num_iterations

# placeholders
mean_vals <- numeric(nr_timepoints * nr_regimes)
lci_vals  <- numeric(nr_timepoints * nr_regimes)
uci_vals  <- numeric(nr_timepoints * nr_regimes)

# Split boot_results_plot into a list of 9 matrices (one per regime)
# Each has 300 (bootstraps) × 29 columns (time points)
reordered_list <- lapply(1:nr_regimes, function(r) {
  boot_results_plot[, seq(from = r, to = ncol(boot_results_plot), by = nr_regimes)]
})

# Compute summary stats 
summary_list <- lapply(reordered_list, function(m) {
  data.frame(
    mean = apply(m, 2, mean),
    lci  = apply(m, 2, quantile, probs = 0.025),
    uci  = apply(m, 2, quantile, probs = 0.975)
  )
})

# Combine results into matrices (rows = time points, cols = regimes)
mean_mat <- sapply(summary_list, \(x) x$mean)
lci_mat  <- sapply(summary_list, \(x) x$lci)
uci_mat  <- sapply(summary_list, \(x) x$uci)

# Turn into data frames
mean_df <- as.data.frame(mean_mat)
lci_df  <- as.data.frame(lci_mat)
uci_df  <- as.data.frame(uci_mat)

# Add row names (T0 to T28)
rownames(mean_df) <- paste0("T", 0:(nr_timepoints-1))
rownames(lci_df)  <- paste0("T", 0:(nr_timepoints-1))
rownames(uci_df)  <- paste0("T", 0:(nr_timepoints-1))

# Add column names (Regime 9 -> 1)
colnames(mean_df) <- paste0("Regime_", 1:nr_regimes)
colnames(lci_df)  <- paste0("Regime_", 1:nr_regimes)
colnames(uci_df)  <- paste0("Regime_", 1:nr_regimes)

stepify <- function(x, y) {
  # Duplicate all but the last x, repeat values to create steps
  x_step <- rep(x, each = 2)[-1]
  y_step <- rep(y, each = 2)[-length(y)*2]
  list(x = x_step, y = y_step)
}

################################################################################

# Open PDF device
pdf(file.path(output_folder_name, "cuminc_suc_ex.pdf"), width = 8, height = 6)  # Landscape

par(mar = c(5, 5, 4, 2) + 0.1)

# Setup base plot
plot(NULL, type = "n", xlim = c(0, 672), ylim = c(0, 1),
     xlab = "Days since time zero", ylab = "Successful extubation",
     xaxt = 'n', yaxt = 'n', cex.lab = 1.5, bty = "n")

# Add horizontal grid lines
abline(h = seq(0.2, 1, by = 0.2), col = "lightgray", lty = 1)

# Axis labels
axis(1, at = c(0, 7*24, 14*24, 21*24, 28*24), labels = c(0, 7, 14, 21, 28), cex.axis = 1.2)
axis(2, at = seq(0, 1, 0.2), las = 1, cex.axis = 1.2)

# prepare time points 
time_points <- seq(0, 672, by = 24)  # 29 time points

# colors for the 3 regimes
colors <- c("#FF9675", "#C04F4C", "#441E1C")

# plot each regime where XlowPEEP == XhighPEEP (regimes 1, 5, 9)
regimes_to_plot <- c(1, 5, 9)

for (idx in seq_along(regimes_to_plot)) {
  
  i <- regimes_to_plot[idx]   # actual regime number
  col <- colors[idx]        # corresponding color
  
  # extract mean, lci, uci for regime i
  mean_vals <- mean_df[, i]
  lci_vals  <- lci_df[, i]
  uci_vals  <- uci_df[, i]
  
  # Stepify all series
  mean_step <- stepify(time_points, mean_vals)
  lci_step  <- stepify(time_points, lci_vals)
  uci_step  <- stepify(time_points, uci_vals)
  
  # Shaded CI (step polygon)
  polygon(
    x = c(lci_step$x, rev(uci_step$x)),
    y = c(lci_step$y, rev(uci_step$y)),
    col = adjustcolor(col, alpha.f = 0.2),
    border = NA
  )
  
  # Mean line as step
  lines(mean_step$x, mean_step$y, col = col, lwd = 2.5, type = "s")
}

legend("bottomright",
       legend = c("P/F>150", "P/F>200", "P/F>250"),
       col = colors, lwd = 2.5, bty = "n",cex = 1.2)

dev.off()
