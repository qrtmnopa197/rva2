// A version of the preregistered model that allows for B_V and alpha to correlate by drawing them from a 
// multivariate normal prior
data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=2> n_f; // Number of fractals (i.e., cues)
  
  array[n_s,n_t] int frac; // The identity of the fractal presented on each trial
  array[n_s,n_t] int alt_frac; // The identity of the alternative fractal on each trial
  
  real out_size; // Absolute value of the binary outcomes

  array[n_s] vector[n_t] rew; // The reward on each trial
  
  int n_vrat; // Total number of valence ratings made
  array[n_s,n_t] int val_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_vrat] val_rat; // The valence rating on each trial
  
  int n_prat; // Total number of probability ratings made
  array[n_s,n_t] int prob_rat_num; // Rating number for each trial
  vector<lower=0,upper=1>[n_prat] prob_rat; // The probability rating on each trial
  
  array[n_s] vector[n_t] prev_vrat_cent; // The previous valence rating on each trial, mean-centered
                                        // Used for autoregressive term
}


parameters{
  // Correlated parameters
  vector[2] bv_alpha_mu; //mean of B_V and of alpha
  vector<lower=0>[2] bv_alpha_sigma; // sd of B_V and alpha
  cholesky_factor_corr[2] bv_alpha_lOmega; // the cholesky factorziation of the correlation matrix between the two
  matrix[2,n_s] bv_alpha_z; // Z-scores representing the deviation of each P's B_V and alpha value from its predicted value

  // Forgetting rate
  real forget_mu;
  real<lower=0> forget_sigma;
  vector[n_s] forget_z;
  
  //Valence predictors
  
  // Intercept
  real B_0_mu;
  real<lower=0> B_0_sigma;
  vector[n_s] B_0_z;
  
  // Effect of reward
  real B_rew_mu; 
  real<lower=0> B_rew_sigma;
  vector[n_s] B_rew_z;
  
  // Effect of alternative cue's V
  real B_alt_mu; 
  real<lower=0> B_alt_sigma;
  vector[n_s] B_alt_z;

  // Autoregressive term
  real B_auto_mu; 
  real<lower=0> B_auto_sigma;
  vector[n_s] B_auto_z;
  
  // Residual SDs
  real<lower=0> prob_sigma;
  
  // Mean valence sigma is implicitly 0 for half-normal prior
  real<lower=0> val_sigma_sigma;
  vector<lower=0>[n_s] val_sigma_z;
  
}

transformed parameters {
  vector[n_prat] prob_pred; // Predicted probability rating
  vector[n_vrat] val_pred; // Predicted valence rating
  vector[n_vrat] val_sigma_ss; // Subject-specific error SDs for valence ratings
  array[n_s,n_t+1] vector[n_f] V; // EVs of each fractals
  {//anonymous scope start
  
    // Get B_V and alpha values for each P
    matrix[2,n_s] bv_alpha = rep_matrix(bv_alpha_mu,n_s) + diag_pre_multiply(bv_alpha_sigma, bv_alpha_lOmega) * bv_alpha_z; 
    // Turn these into vectors of subject-level parameters, as you more typically work with
    vector[n_s] B_V = bv_alpha[1,]';
    vector[n_s] alpha = inv_logit(bv_alpha[2,]');
  
    // Get subject-specific values based on non-centered parameterization
    vector[n_s] forget = inv_logit(forget_mu + forget_sigma*forget_z);
    
    vector[n_s] val_sigma = val_sigma_sigma*val_sigma_z;
    
    vector[n_s] B_0 = B_0_mu + B_0_sigma*B_0_z;
    vector[n_s] B_rew = B_rew_mu + B_rew_sigma*B_rew_z;
    vector[n_s] B_alt = B_alt_mu + B_alt_sigma*B_alt_z;
    vector[n_s] B_auto = B_auto_mu + B_auto_sigma*B_auto_z;
    
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          V[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all V-values to 0
        }
        V[s,t+1] = V[s,t]*forget[s]; // Decay toward 0 by defafult
        V[s,t+1,frac[s,t]] = V[s,t,frac[s,t]] + alpha[s]*(rew[s,t] - V[s,t,frac[s,t]]); // Update V of presented cue
        
        // If the participant made a probability rating...
        if(prob_rat_num[s,t] != 0){
          // Add to the vector of predicted Vs - first transforming V to be on the probability scale
          prob_pred[prob_rat_num[s,t]] = 0.5 + V[s,t+1,frac[s,t]]/(2*out_size); 
        }
        // If the participant made a valence rating... 
        if(val_rat_num[s,t] != 0){
          // Add to the vector of predicted ratings
          val_pred[val_rat_num[s,t]] = inv_logit(B_0[s] + B_rew[s]*rew[s,t] + B_V[s]*V[s,t,frac[s,t]] 
                                                + B_alt[s]*V[s,t,alt_frac[s,t]] + B_auto[s]*prev_vrat_cent[s,t]); 
          val_sigma_ss[val_rat_num[s,t]] = val_sigma[s]; //the subject-specific error SD for this observation
        }
      }
    }
  }//anonymous scope end
}
model{
  // Set priors on B_V and alpha means and sds, along with the correlation between them
  bv_alpha_mu[1] ~ normal(0,0.02); 
  bv_alpha_sigma[1] ~ normal(0,0.03);
  bv_alpha_mu[2] ~ normal(-.05,1.7);
  bv_alpha_sigma[2] ~ normal(0,4);
  to_vector(bv_alpha_z) ~ std_normal();
  bv_alpha_lOmega ~ lkj_corr_cholesky(1);
  
  forget_mu ~ normal(-.05,1.7); //approx. equal to uniform distribution when passed through inv_logit
  forget_sigma ~ normal(0,4);
  forget_z ~ std_normal();
  
  B_0_mu ~ normal(0,2);
  B_0_sigma ~ normal(0,3);
  B_0_z ~ std_normal();

  B_rew_mu ~ normal(0,0.03); 
  B_rew_sigma ~ normal(0,0.05);
  B_rew_z ~ std_normal();
  
  B_alt_mu ~ normal(0,0.03); 
  B_alt_sigma ~ normal(0,0.05);
  B_alt_z ~ std_normal();

  B_auto_mu ~ normal(0,2); 
  B_auto_sigma ~ normal(0,3);
  B_auto_z ~ std_normal();
  
  prob_sigma ~ normal(0,1);
  
  val_sigma_sigma ~ normal(0,2);
  val_sigma_z ~ std_normal();
  
  val_rat ~ normal(val_pred, val_sigma_ss) T[0,1];
  prob_rat ~ normal(prob_pred, prob_sigma) T[0,1];
}
generated quantities{
  real bv_alpha_cor = bv_alpha_lOmega[2,1]; // correlation between B_V and alpha
  vector[n_vrat] val_lik; //log likelihoods of all valence ratings
  vector[n_prat] prob_lik; //log likelihoods of all probability ratings
  
  for(v in 1:n_vrat){
    val_lik[v] = normal_lpdf(val_rat[v] | val_pred[v], val_sigma_ss[v]);
  }
  for(p in 1:n_prat){
    prob_lik[p] = normal_lpdf(prob_rat[p] | prob_pred[p], prob_sigma);
  }
}
