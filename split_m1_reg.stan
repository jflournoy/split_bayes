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
  initV  = rep_vector(0.0, K);
}
parameters { # Need to add cue-dimension to mu and individual priors
  # declare as vectors for vectorizing
  real mu_xi;  
  vector[D] mu_ep;
  real mu_rho;
  real<lower=0> sigma_xi;
  corr_matrix[D] Omega_ep; // prior correlation
  vector<lower=0>[D] tau_ep; // prior scale
  real<lower=0> sigma_b;
  real<lower=0> sigma_rho;
  vector[N] xi_pr;         # noise 
  vector[D] ep_pr[N];      # learning rate 
  vector[N] rho_pr;        # rho, inv temp 
} 
transformed parameters{
  vector<lower=0,upper=1>[N] xi;
  vector<lower=0,upper=1>[D] ep[N];
  vector<lower=0>[N] rho;
  cov_matrix[D] sigma_ep;
  matrix[D, D] L_ep;
  
  sigma_ep =  quad_form_diag(Omega_ep, tau_ep);
  L_ep = cholesky_decompose(sigma_ep);
  
  for (i in 1:N) {
    xi[i]  = Phi_approx( mu_xi + sigma_xi * xi_pr[i] );
    ep[i]  = Phi_approx( mu_ep + L_ep * ep_pr[i] );
  }
  rho = exp( mu_rho + sigma_rho * rho_pr );
}
model {  
# gng_m2: RW + noise + bias model in Guitart-Masip et al 2012
  # hyper parameters
  mu_xi  ~ normal(0, 1.0); 
  mu_ep  ~ normal(0, 1.0);
  mu_rho  ~ normal(0, 1.0); 
  
  sigma_xi ~ cauchy(0, 5.0);
  tau_ep ~ cauchy(0, 2.5);
  Omega_ep ~ lkj_corr(2);
  sigma_rho ~ cauchy(0, 5.0);
  
  # individual parameters w/ Matt trick
  xi_pr  ~ normal(0, 1.0);
  for (i in 1:N) {
    ep_pr[i] ~ normal(0, 1.0);
  }
  rho_pr ~ normal(0, 1.0);

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
      pR[ cue[i,t] ]   = pR[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
      pressed_r[i,t] ~ bernoulli( pR[ cue[i,t] ] );
      
      # update action values
      if (pressed_r[i,t]) { # update go value 
        qv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[i] * outcome[i,t] - qv_r[ cue[i,t] ]);
      } else { # update no-go value  
        qv_l[ cue[i,t] ] = qv_l[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[i] * outcome[i,t] - qv_l[ cue[i,t] ]);  
      }  
    } # end of t loop
  } # end of i loop
}
generated quantities {
  real<lower=0, upper=1> mu_xi_o;
  vector<lower=0, upper=1>[D] mu_ep_o;
  real<lower=0> mu_rho_o;
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
        pR[ cue[i,t] ]   = pR[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
        log_lik[i] = log_lik[i] + bernoulli_lpmf( pressed_r[i,t] | pR[ cue[i,t] ] );
        
        # Model regressors --> store values before being updated
        Qright[i,t]   = qv_r[ cue[i,t] ];
        Qleft[i,t] = qv_l[ cue[i,t] ];
        Wright[i,t]   = wv_r[ cue[i,t] ];
        Wleft[i,t] = wv_l[ cue[i,t] ];
        
        # update action values
        if (pressed_r[i,t]) { # update go value 
          qv_r[ cue[i,t] ]  = qv_r[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[i] * outcome[i,t] - qv_r[ cue[i,t] ]);
        } else { # update no-go value  
          qv_l[ cue[i,t] ] = qv_l[ cue[i,t] ] + ep[ i, cuetype[i,t] ] * (rho[i] * outcome[i,t] - qv_l[ cue[i,t] ]);  
        }  
      } # end of t loop
    } # end of i loop
  } # end of local section
}
