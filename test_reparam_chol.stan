# _reg: generates model-based regressors
data {
  int<lower=1> N; #Number of subjects
  int<lower=1> T; #Max number of trials
  int<lower=1> K; #Number of distinct cues
  int<lower=1> D; #Number of groups of cues
  int<lower=1, upper=T> Tsubj[N]; #number of trials per subject
  int<lower=0, upper=1> pressed_r[N, T]; #Subject x Trial matrix of button presses (L = 0, R = 1)
  int<lower=1, upper=K> cue[N, T]; #Subject x Trial matrix of which cue was presented
  int<lower=1, upper=D> cuetype[N, T]; #Subject x Trial matrix of which type of cue was presented
}
transformed data{
  matrix[N,1] u;
  u[,1] = rep_vector(1,N);
}
parameters { # Need to add cue-dimension to mu and individual priors
  matrix[D, N] ep_pr; //individual deviations "z" p149
  cholesky_factor_corr[D] L_Omega_ep; // prior correlation
  vector<lower=0,upper=pi()/2>[D] tau_unif_ep; // prior scale
  matrix[1,D] mu_ep; //group-level intercept "gamma" p149
} 
transformed parameters{
  matrix<lower=0,upper=1>[N,D] ep;
  vector<lower=0>[D] tau_ep; //prior scale p 149

  for (d in 1:D) tau_ep[d] = 2.5 * tan(tau_unif_ep[d]);
  
  ep  = Phi_approx( u*mu_ep + (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr )' );
}
model {  
  to_vector(mu_ep) ~ normal(0, 1); 
  L_Omega_ep ~ lkj_corr_cholesky(2);
  to_vector(ep_pr) ~ normal(0, .5);
  # individual parameters w/ Matt trick

  for (i in 1:N) {
    for (t in 1:Tsubj[i])  {
      vector[K] pR;   # prob of R (press)
      pR[ cue[i,t] ] = inv_logit( ep[ i, cuetype[i,t] ] );
      pressed_r[i,t] ~ bernoulli( pR[ cue[i,t] ] );
    } # end of t loop
  } # end of i loop
}
generated quantities {
  matrix[N,D] ep_raw;
  matrix[D,D] some_sig;
  matrix[N,D] indiv_var;
  matrix[D,D] other_sig;
  
  indiv_var = (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr )';
  ep_raw = u*mu_ep + (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr )';
  some_sig = diag_pre_multiply(tau_ep,L_Omega_ep);
  other_sig = (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr ) * (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr )';
}
