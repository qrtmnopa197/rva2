power_rva_2 <- function(blocks = 2,
                        trials = 40,
                        min_prob = .25,
                        max_prob =.75,
                        drift = 0.05,
                        loss = -100,
                        win = 100,
                        sims = 500,
                        sigma = .12,
                        sample_sizes = c(40,50,60,70)){
  
  return_df <- data.frame("sample_size"=c(),"power_reward"=c(),"power_rpe"=c(),"reward_beta"=c(),"rpe_beta"=c(),"rsq"=c())
  for(ss in sample_sizes){
    rew_vec <- c() #number of times the effect for reward has been significant
    rew_betas <- c() #standardized effect of reward on valence
    rpe_vec <- c() #ditto RPE
    rpe_betas <- c()
    rsq <- c() # the R-squared value for the regression model on each simulation
    for(i in 1:sims){
      #initialize df containing data for this simulation
      sim_data_df <- data.frame("val"=c(),"rew"=c(),"rpe"=c(),"avv"=c(),"sub"=c()) 
      for(s in 1:ss){
        #get parameters for this subject
        alpha <- sample(c(.10,.15,.20,.25,.30),prob=c(.05,.25,.4,.25,.05),size=1)
        forget <- sample(c(.10,.15,.20,.25,.30),prob=c(.05,.25,.4,.25,.05),size=1)
        b_rew <- rnorm(1,mean=0.007,sd=.0035)
        b_rpe <- rnorm(1,mean=0.0012,sd=.0006)
        b_avv <- rnorm(1,mean=0.005,sd=0.0025)
        b_pr <- rnorm(1,mean=0.4,sd=0.2)
        
        for(b in 1:blocks){
          #initialize vectors that will go in the sim_data_df
          val <- rep(NA,trials)
          rew <- rep(0,trials)
          rpe <- rep(0,trials)
          avv <- rep(0,trials)
          pr <- rep(0,trials)
          prev_rat <- 0
          
          # Initialize V
          V <- c(0,0)
          
          # Determine the trial type on each trial
          # If it's less than 20, it's cue 1
          # If it's even, it's an affect rating trial
          tr_typ <- c()
          for(i in 1:(trials/4)){
            tr_typ <- c(tr_typ,sample(c(10,11,20,21)))
          }
          
          
          # Initialize likelihoods of winning from each cue
          win_prob <- c(runif(1, min_prob, max_prob),runif(1, min_prob, max_prob))
          
          for(t in 1:trials){
            # Determine whether cue 1 or 2 was presented
            if(tr_typ[t] < 20){
              c <- 1
              u <- 2
            } else{
              c <- 2
              u <- 1
            }
            # Generate outcome
            rand <- runif(1)
            if(rand < win_prob[c]){
              rew[t] <- win
            } else{
              rew[t] <- loss
            }
            # Calculate valence predictors
            rpe[t] <- rew[t] - V[c]
            avv[t] <- sum(V)/2
            pr[t] <- prev_rat - 0.5
            if(tr_typ[t] %% 2 == 0){
              # Simulate valence rating
              val[t] <- rnorm(1,sigmoid(b_rew*rew[t] + b_rpe*rpe[t] 
                                        + b_avv*avv[t] + b_pr*pr[t]), sigma) 
              prev_rat <- val[t]
            }
            
            
            V[c] <- V[c] + alpha*rpe[t] # Update V (for next trial) based on the current trial's outcome
            V[u] <- V[u]*(1-forget)
            
            #drift probability of winning
            win_prob[c] <- win_prob[c] + rnorm(1,0,drift)
            if(win_prob[c] < min_prob){
              win_prob[c] <- min_prob + (min_prob - win_prob[c])
            } else{
              if(win_prob[c] > max_prob){
                win_prob[c] <- max_prob - (win_prob[c] - max_prob)
              }
            }
          }
          # Add this block's data to the general DF for this simulation
          block_df <- data.frame("val" = val,"rew" = rew, "rpe" = rpe, "avv" = avv, "pr" = pr, "sub" = s) 
          sim_data_df <- rbind(sim_data_df,block_df)
        }
      }
      #Fit MLM to the data. If there is an error, rerun the MLM without random effects, as errors typically occur 
      #when the random effects structure is overly complicated and there are convergence issues.
      sim_fit <- tryCatch(
        {lme(val ~ rew + rpe + avv + pr,random =~ 0 + rew + rpe + avv + pr | sub, sim_data_df, na.action=na.omit)},
        error = function(e) {lm(val ~ rew + rpe + avv + pr, sim_data_df, na.action=na.omit)}
      )
      sum_sim_fit <- summary(sim_fit)
      rat_sdd <- filter(sim_data_df,!is.na(val))
      #Record whether the reward and RPE coefficients were significant, and record effect sizes
      if(class(sim_fit) == "lm"){
        B_rew <- sum_sim_fit$coefficients[2,1]*(sd(rat_sdd$rew)/sd(rat_sdd$val))
        B_rpe <- sum_sim_fit$coefficients[3,1]*(sd(rat_sdd$rpe)/sd(rat_sdd$val))
        rew_betas <- c(rew_betas,B_rew)
        rpe_betas <- c(rpe_betas,B_rpe)
        rsq <- c(rsq,sum_sim_fit$r.squared)
        
        if(sum_sim_fit$coefficients[2,4] < .10 && B_rew > 0){
          rew_vec <- c(rew_vec,1)
        } else{
          rew_vec <- c(rew_vec,0)
        }
        if(sum_sim_fit$coefficients[3,4] < .10 && B_rpe > 0){
          rpe_vec <- c(rpe_vec,1)
        } else{
          rpe_vec <- c(rpe_vec,0)
        }
      } else if(class(sim_fit) == "lme"){
        B_rew <- sum_sim_fit$tTable[2,1]*(sd(rat_sdd$rew)/sd(rat_sdd$val))
        B_rpe <- sum_sim_fit$tTable[3,1]*(sd(rat_sdd$rpe)/sd(rat_sdd$val))
        rew_betas <- c(rew_betas,B_rew)
        rpe_betas <- c(rpe_betas,B_rpe)
        rsq <- c(rsq,r.squaredGLMM(sim_fit))
        
        if(sum_sim_fit$tTable[2,5] < .10 && B_rew > 0){
          rew_vec <- c(rew_vec,1)
        } else{
          rew_vec <- c(rew_vec,0)
        }
        if(sum_sim_fit$tTable[3,5] < .10 && B_rpe > 0){
          rpe_vec <- c(rpe_vec,1)
        } else{
          rpe_vec <- c(rpe_vec,0)
        }
      }
    }
    return_df <- rbind(return_df,
                       data.frame("sample_size"=ss,"power_reward"=mean(rew_vec),"power_rpe"=mean(rpe_vec),
                                  "reward_beta"=mean(rew_betas),"rpe_beta"=mean(rpe_betas),"rsq"=mean(rsq))
    )
  }
  return(return_df)
}


