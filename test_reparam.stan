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
transformed data {
  vector[K] initV;
  initV  = rep_vector(0.0, K);
}
parameters { # Need to add cue-dimension to mu and individual priors
  # declare as vectors for vectorizing
  vector[D] mu_ep;
  corr_matrix[D] Omega_ep; // prior correlation
  vector<lower=0>[D] tau_ep; // prior scale
  vector[D] ep_pr[N];
} 
transformed parameters{
  vector<lower=0,upper=1>[D] ep[N];
  cov_matrix[D] sigma_ep;
  matrix[D, D] L_ep;
  
  sigma_ep =  quad_form_diag(Omega_ep, tau_ep);
  L_ep = cholesky_decompose(sigma_ep);
  
  for (i in 1:N) {
    ep[i]  = Phi_approx( mu_ep + L_ep * ep_pr[i] );
  }
}
model {  
  # gng_m2: RW + noise + bias model in Guitart-Masip et al 2012
  # hyper parameters

  mu_ep  ~ normal(0, 1.0); 

  tau_ep ~ cauchy(0, 2.5);
  Omega_ep ~ lkj_corr(2);

  # individual parameters w/ Matt trick
  for (i in 1:N) {
    ep_pr[i] ~ normal(0, 1.0);
  }
  
  for (i in 1:N) {
    for (t in 1:Tsubj[i])  {
      vector[K] pR;   # prob of R (press) 
      pR[ cue[i,t] ]   = inv_logit( ep[ i, cuetype[i,t] ] ); 
      pressed_r[i,t] ~ bernoulli( pR[ cue[i,t] ] );
    } # end of t loop
  } # end of i loop
}
generated quantities {
  cov_matrix[D] sigma_ep_o;
  matrix[D, D] L_ep_o;
  vector<lower=0, upper=1>[D] mu_ep_o;
  
  sigma_ep_o =  quad_form_diag(Omega_ep, tau_ep);
  L_ep_o = cholesky_decompose(sigma_ep_o);
  mu_ep_o  = Phi_approx(mu_ep);
}