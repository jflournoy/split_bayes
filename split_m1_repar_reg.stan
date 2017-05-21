# _reg: generates model-based regressors
data {
  int<lower=1> N; #Number of subjects
  int<lower=1> T; #Max number of trials
  int<lower=1> K; #Number of distinct cues
  int<lower=1> D; #Number of groups of cues
  int<lower=1, upper=T> Tsubj[N]; #number of trials per subject
  real outcome[N, T]; #Subject x Trial matrix of outcomes
  int<lower=0, upper=1> pressed_r[N, T]; #Subject x Trial matrix of button presses (L = 0, R = 1)
  int<lower=1, upper=K> cue[N, T]; #Subject x Trial matrix of which cue was presented
  int<lower=1, upper=D> cuetype[N, T]; #Subject x Trial matrix of which type of cue was presented
}
transformed data {
  vector[K] initV;  
  matrix[N,1] u;
  
  initV  = rep_vector(0.0, K);
  u[,1] = rep_vector(1,N);
}
parameters { # Need to add cue-dimension to mu and individual priors
  # declare as vectors for vectorizing
  
  #Learning rate
  matrix[D, N] ep_pr; //individual deviations "z" p149
  cholesky_factor_corr[D] L_Omega_ep; // prior correlation
  vector<lower=0,upper=pi()/2>[D] tau_unif_ep; // prior scale
  matrix[1,D] mu_ep; //group-level intercept "gamma" p149
  
  #Noise
  matrix[D, N] xi_pr; //individual deviations "z" p149
  cholesky_factor_corr[D] L_Omega_xi; // prior correlation
  vector<lower=0,upper=pi()/2>[D] tau_unif_xi; // prior scale
  matrix[1,D] mu_xi; //group-level intercept "gamma" p149
  
  #Inverse Temp
  matrix[D, N] rho_pr; //individual deviations "z" p149
  cholesky_factor_corr[D] L_Omega_rho; // prior correlation
  vector<lower=0,upper=pi()/2>[D] tau_unif_rho; // prior scale
  matrix[1,D] mu_rho; //group-level intercept "gamma" p149 
} 
transformed parameters{
  matrix<lower=0,upper=1>[N,D] ep;
  vector<lower=0>[D] tau_ep; //prior scale p 149

  matrix<lower=0,upper=1>[N,D] xi;
  vector<lower=0>[D] tau_xi; //prior scale p 149
  
  matrix<lower=0,upper=1>[N,D] rho;
  vector<lower=0>[D] tau_rho; //prior scale p 149

  for (d in 1:D) tau_ep[d] = 2.5 * tan(tau_unif_ep[d]);
  for (d in 1:D) tau_xi[d] = 2.5 * tan(tau_unif_xi[d]);
  for (d in 1:D) tau_rho[d] = 2.5 * tan(tau_unif_rho[d]);

  xi  = Phi_approx( u*mu_xi + (diag_pre_multiply(tau_xi,L_Omega_xi) * xi_pr )' );
  ep  = Phi_approx( u*mu_ep + (diag_pre_multiply(tau_ep,L_Omega_ep) * ep_pr )' );
  rho = exp( u*mu_rho + (diag_pre_multiply(tau_rho,L_Omega_rho) * rho_pr )' );
}
model {  
  # gng_m2: RW + noise + bias model in Guitart-Masip et al 2012
  # hyper parameters
  to_vector(mu_ep) ~ normal(0, 1); 
  L_Omega_ep ~ lkj_corr_cholesky(2);
  to_vector(ep_pr) ~ normal(0, .5);
  
  to_vector(mu_xi) ~ normal(0, 1); 
  L_Omega_xi ~ lkj_corr_cholesky(2);
  to_vector(xi_pr) ~ normal(0, 1);
  
  to_vector(mu_rho) ~ normal(0, 1); 
  L_Omega_rho ~ lkj_corr_cholesky(2);
  to_vector(rho_pr) ~ normal(0, 1);

  for (i in 1:N) {
    vector[K] wv_r;  # action weight for r
    vector[K] wv_l; # action weight for l
    vector[K] qv_r;  # Q value for r
    vector[K] qv_l; # Q value for l
    vector[K] pR;   # prob of R (press) 

    wv_r  = initV;
    wv_l = initV;
    qv_r  = initV;
    qv_l = initV;
  
    for (t in 1:Tsubj[i])  {
      wv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ];
      wv_l[ cue[i,t] ] = qv_l[ cue[i,t] ];  # qv_l is always equal to wv_r (regardless of action)      
      pR[ cue[i,t] ]   = inv_logit( wv_r[ cue[i,t] ] - wv_l[ cue[i,t] ] ); 
      pR[ cue[i,t] ]   = pR[ cue[i,t] ] * (1 - xi[ i, cuetype[i,t] ]) + xi[ i, cuetype[i,t] ]/2;  # noise
      pressed_r[i,t] ~ bernoulli( pR[ cue[i,t] ] );
      
      # update action values
      if (pressed_r[i,t]) { # update go value 
        qv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[ i, cuetype[i,t] ] * outcome[i,t] - qv_r[ cue[i,t] ]);
      } else { # update no-go value  
        qv_l[ cue[i,t] ] = qv_l[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[ i, cuetype[i,t] ] * outcome[i,t] - qv_l[ cue[i,t] ]);  
      }  
    } # end of t loop
  } # end of i loop
}
generated quantities {
  matrix<lower=0, upper=1>[1,D] mu_xi_o;
  matrix<lower=0, upper=1>[1,D] mu_ep_o;
  matrix<lower=0>[1,D] mu_rho_o;
  real log_lik[N];
  real Qright[N, T];
  real Qleft[N, T];
  real Wright[N, T];
  real Wleft[N, T];
  
  mu_xi_o  = Phi_approx(mu_xi);
  mu_ep_o  = Phi_approx(mu_ep);
  mu_rho_o = exp(mu_rho); 
  
  { # local section, this saves time and space
    for (i in 1:N) {
      vector[K] wv_r;  # action wegith for go
      vector[K] wv_l; # action wegith for nogo
      vector[K] qv_r;  # Q value for go
      vector[K] qv_l; # Q value for nogo
      vector[K] pR;   # prob of go (press) 
  
      wv_r  = initV;
      wv_l = initV;
      qv_r  = initV;
      qv_l = initV;
    
      log_lik[i] = 0;

      for (t in 1:T)  {
        wv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ];
        wv_l[ cue[i,t] ] = qv_l[ cue[i,t] ];  # qv_ng is always equal to wv_ng (regardless of action)      
        pR[ cue[i,t] ]   = inv_logit( wv_r[ cue[i,t] ] - wv_l[ cue[i,t] ] ); 
        pR[ cue[i,t] ]   = pR[ cue[i,t] ] * (1 - xi[ i, cuetype[i,t] ]) + xi[ i, cuetype[i,t] ]/2;  # noise
        log_lik[i] = log_lik[i] + bernoulli_lpmf( pressed_r[i,t] | pR[ cue[i,t] ] );
        
        # Model regressors --> store values before being updated
        Qright[i,t]   = qv_r[ cue[i,t] ];
        Qleft[i,t] = qv_l[ cue[i,t] ];
        Wright[i,t]   = wv_r[ cue[i,t] ];
        Wleft[i,t] = wv_l[ cue[i,t] ];
        
        # update action values
        if (pressed_r[i,t]) { # update go value 
          qv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[ i, cuetype[i,t] ] * outcome[i,t] - qv_r[ cue[i,t] ]);
        } else { # update no-go value  
          qv_l[ cue[i,t] ] = qv_l[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[ i, cuetype[i,t] ] * outcome[i,t] - qv_l[ cue[i,t] ]);  
        }  
      } # end of t loop
    } # end of i loop
  } # end of local section
}
