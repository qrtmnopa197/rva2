# Function to get the nth previous outcome for a given frac_ix
get_lagged_outcome <- function(data, match_from_col, match_to_col, lag_n) {
  sapply(seq_len(nrow(data)), function(i) {
    if (i == 1) return(0)
    match_value <- data[[match_from_col]][i]
    prior_rows <- which(data[[match_to_col]][1:(i - 1)] == match_value)
    if (length(prior_rows) >= lag_n) {
      data$outcome[prior_rows][length(prior_rows) - lag_n + 1]
    } else {
      0
    }
  })
}

#Returns a list of data for input to stan, given trial-level data
#trials: trial-level data
#n_t: number of trials; must be set manually
stan_data_rva2 <- function(trials,n_t){
  n_s <- length(unique(trials$id)) #get number of subjects
  n_f <- max(trials$frac_ix)
  
  out_size <- abs(trials$outcome[1]) # Absolute value of outcomes
  
  rew <- sub_by_trial_vec_list(trials,"outcome") # Outcome on each trial
  
  frac <- sub_by_trial_matrix(trials,"frac_ix") # Fractal index on each trial
  alt_frac <- sub_by_trial_matrix(trials,"alt_frac_ix") # Alternative fractal on each trial
  
  #valence rating data
  val_rat_num <- sub_by_trial_matrix(trials,"vrat_number")
  n_vrat <- max(trials$vrat_number)
  
  val_rat_trials <- filter(trials,vrat_number != 0)
  val_rat <- val_rat_trials$val_rat
  
  #probability rating data
  prob_rat_num <- sub_by_trial_matrix(trials,"prat_number")
  n_prat <- max(trials$prat_number)
  
  prob_rat_trials <- filter(trials,prat_number != 0)
  prob_rat <- prob_rat_trials$prob_rat
  
  prev_vrat_cent <- sub_by_trial_vec_list(trials,"prev_vrat_cent") #previous valence rating
  
  data <- list(
    n_t = n_t,
    n_s = n_s,
    n_f = n_f,
    out_size = out_size,
    rew = rew,
    frac = frac,
    alt_frac = alt_frac,
    val_rat_num = val_rat_num,
    n_vrat = n_vrat,
    val_rat = val_rat,
    prob_rat_num = prob_rat_num,
    n_prat = n_prat,
    prob_rat = prob_rat,
    prev_vrat_cent = prev_vrat_cent
  )
  return(data)
}


# Add multiple columns to a df with the lagged outcome of the currently shown and currently unshown fractal
# block_df: block-level data frame
# lag: how many trials back to lag
# out_col: name of the column with the outcomes
add_mult_shunsh_lag_outs <- function(block_df,lags,out_col = "outcome"){
  for(lag in lags){
    block_df <- add_shunsh_lag_out(block_df,lag,out_col=out_col)
  }
  return(block_df)
}

# Adds  columns to a df with the lagged outcome of the currently shown and currently unshown fractal
# block_df: block-level data frame
# lag: how many trials back to lag
# out_col: name of the column with the outcomes
add_shunsh_lag_out <- function(block_df,lag,out_col = "outcome"){
  # Create columns for the outcome of the un/shown cue, lagged the appropriate number of trials
  sh_col <- paste0("out_shc_lag",lag)
  unsh_col <- paste0("out_unshc_lag",lag)
  block_df[c(sh_col,unsh_col)] <- 0
  
  # For each trial...
  for(t in (lag+1):nrow(block_df)){
    if(block_df$frac_ix[t] == block_df$frac_ix[t-lag]){
      block_df[t,sh_col] <- block_df[t-lag,out_col]
    } else{
      block_df[t,unsh_col] <- block_df[t-lag,out_col]
    }
  }
  return(block_df)
}