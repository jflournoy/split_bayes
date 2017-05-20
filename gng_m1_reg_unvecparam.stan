# _reg: generates model-based regressors
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=2> K;
  int<lower=1, upper=T> Tsubj[N];
  real outcome[N, T];
  int<lower=0, upper=1> pressed[N, T];
  int<lower=1, upper=K> cue[N, T];
}
transformed data {
  vector[K] initV;
  initV  = rep_vector(0.0, K);
}
parameters {
  # declare as vectors for vectorizing
  real mu_xi;  
  real mu_ep;
  real mu_rho;  
  real<lower=0> sigma_xi;
  real<lower=0> sigma_ep;
  real<lower=0> sigma_rho;
  vector[N] xi_pr;          # noise 
  vector[N] ep_pr;          # learning rate 
  vector[N] rho_pr;         # rho, inv temp 
}
transformed parameters{
  vector<lower=0,upper=1>[N] xi;
  vector<lower=0,upper=1>[N] ep;
  vector<lower=0>[N] rho;
     
  for (i in 1:N) {
    xi[i]  = Phi_approx( mu_xi + sigma_xi * xi_pr[i] );
    ep[i]  = Phi_approx( mu_ep + sigma_ep * ep_pr[i] );
  }
  rho = exp( mu_rho + sigma_rho * rho_pr );
}
model {  
# gng_m1: RW + noise model in Guitart-Masip et al 2012
  # hyper parameters
  mu_xi  ~ normal(0, 1.0); 
  mu_ep  ~ normal(0, 1.0);
  mu_rho  ~ normal(0, 1.0);  
  sigma_xi ~ cauchy(0, 5.0);
  sigma_ep ~ cauchy(0, 5.0);
  sigma_rho ~ cauchy(0, 5.0);
  
  # individual parameters w/ Matt trick
  xi_pr  ~ normal(0, 1.0);   
  ep_pr  ~ normal(0, 1.0);   
  rho_pr ~ normal(0, 1.0);

  for (i in 1:N) {
    vector[K] wv_g;  # action weight for go
    vector[K] wv_ng; # action weight for nogo
    vector[K] qv_g;  # Q value for go
    vector[K] qv_ng; # Q value for nogo
    vector[K] pGo;   # prob of go (press) 

    wv_g  = initV;
    wv_ng = initV;
    qv_g  = initV;
    qv_ng = initV;
  
    for (t in 1:Tsubj[i])  {
      wv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ];
      wv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ];  # qv_ng is always equal to wv_ng (regardless of action)      
      pGo[ cue[i,t] ]   = inv_logit( wv_g[ cue[i,t] ] - wv_ng[ cue[i,t] ] ); 
      pGo[ cue[i,t] ]   = pGo[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
      pressed[i,t] ~ bernoulli( pGo[ cue[i,t] ] );
      
      # update action values
      if (pressed[i,t]) { # update go value 
        qv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + ep[i] * (rho[i] * outcome[i,t] - qv_g[ cue[i,t] ]);
      } else { # update no-go value  
        qv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ] + ep[i] * (rho[i] * outcome[i,t] - qv_ng[ cue[i,t] ]);  
      }  
    } # end of t loop
  } # end of i loop
}
generated quantities {
  real<lower=0, upper=1> mu_xi_o;
  real<lower=0, upper=1> mu_ep_o;
  real<lower=0> mu_rho_o;
  real log_lik[N];
  real Qgo[N, T];
  real Qnogo[N, T];
  real Wgo[N, T];
  real Wnogo[N, T];

  mu_xi_o  = Phi_approx(mu_xi);
  mu_ep_o  = Phi_approx(mu_ep);
  mu_rho_o = exp(mu_rho); 
  
  { # local section, this saves time and space
    for (i in 1:N) {
      vector[K] wv_g;  # action wegith for go
      vector[K] wv_ng; # action wegith for nogo
      vector[K] qv_g;  # Q value for go
      vector[K] qv_ng; # Q value for nogo
      vector[K] pGo;   # prob of go (press) 
  
      wv_g  = initV;
      wv_ng = initV;
      qv_g  = initV;
      qv_ng = initV;
    
      log_lik[i] = 0;

      for (t in 1:T)  {
        wv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ];
        wv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ];  # qv_ng is always equal to wv_ng (regardless of action)      
        pGo[ cue[i,t] ]   = inv_logit( wv_g[ cue[i,t] ] - wv_ng[ cue[i,t] ] ); 
        pGo[ cue[i,t] ]   = pGo[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
        log_lik[i] = log_lik[i] + bernoulli_lpmf( pressed[i,t] | pGo[ cue[i,t] ] );
        
        # Model regressors --> store values before being updated
        Qgo[i,t]   = qv_g[ cue[i,t] ];
        Qnogo[i,t] = qv_ng[ cue[i,t] ];
        Wgo[i,t]   = wv_g[ cue[i,t] ];
        Wnogo[i,t] = wv_ng[ cue[i,t] ];
        
        # update action values
        if (pressed[i,t]) { # update go value 
          qv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + ep[i] * (rho[i] * outcome[i,t] - qv_g[ cue[i,t] ]);
        } else { # update no-go value  
          qv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ] + ep[i] * (rho[i] * outcome[i,t] - qv_ng[ cue[i,t] ]);  
        }  
      } # end of t loop
    } # end of i loop
  } # end of local section
}
