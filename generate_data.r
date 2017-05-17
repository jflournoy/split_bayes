#' # Generate data in Stan
#' 
#' Just a quick example of how one can generate data in stan

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model <- '
data {
  int<lower=1> N;
  real x[N];
}
transformed data {
  vector[N] mu;
  cov_matrix[N] Sigma;
  for (i in 1:N){
    mu[i] = 0;
  }  
  for (i in 1:N)
    for (j in 1:N)
      Sigma[i, j] = exp(-pow(x[i] - x[j], 2)) + (i == j ? 0.1 : 0.0);
  print(Sigma)
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu, Sigma);
}
generated quantities {
  cov_matrix[N] Sigma_out;
  Sigma_out = Sigma;
}'

N=5

somedata <- list(N=N,
                 x=runif(N,0,1))

fit <- stan(model_code = model, data = somedata,  iter = 10000, chains = 8)
fit <- stan(model_code = model, data = somedata,  iter = 1, chains = 1)

list_of_draws <- extract(fit, permuted = F)

plot(fit)
stan_trace(fit)
fit@par_dims

lapply(1:8, function(x){
  cov(list_of_draws[,x,1:5])
})

list_of_draws[1,8,6:(6+24)]

list_of_draws_perm <- extract(fit, permuted = T)

cov(list_of_draws$y)

#' # Generate binomial data
#' 
#' This will assume y generated above is on the logit scale and use it as the predictor matrix to generate a probability for a 
#' bernoulli variable.

bin_model <- '
data {
  int<lower=1> N;
  real x[N];
  vector[N] beta;
}
transformed data {
  vector[N] mu;
  cov_matrix[N] Sigma;
  for (i in 1:N){
    mu[i] = 0;
}  
for (i in 1:N)
  for (j in 1:N)
    Sigma[i, j] = exp(-pow(x[i] - x[j], 2)) + (i == j ? 0.1 : 0.0);
print(Sigma)
}
parameters {
  matrix[1,N] X;
}
model {
  X[1,] ~ multi_normal(mu, Sigma);
}
generated quantities {
  real y_out;
  real plogit;
  plogit = X[1,]*beta;
  y_out = bernoulli_logit_rng(plogit);
}'

N=5

bin_data <- list(N=N,
                 x=runif(N,0,1),
                 beta=seq(.1,.5, length.out = N))

bin_fit <- stan(model_code = bin_model, data = bin_data,  iter = 1000, chains = 8)

bin_fit

#' # Generate data like RW model
#' 
#' 

rw_model <- '
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> K;
  int<lower=1, upper=T> Tsubj[N];
  real outcome_l[N, T];
  real outcome_r[N, T];
  #int<lower=0, upper=1> pressed_r[N, T];
  int<lower=1, upper=K> cue[N, T];
}
transformed data {
  vector[K] initV;
  vector[4] mu_p; 
  initV  = rep_vector(0.0, K);
  mu_p[1]  = 0;
  mu_p[2]  = .5;
  mu_p[3]  = 0; 
  mu_p[4]  = 1; 
}
parameters {
  # declare as vectors for vectorizing
  # vector[4] mu_p;  
  vector<lower=0>[4] sigma; 
  vector[N] xi_pr;         # noise 
  vector[N] ep_pr;         # learning rate 
  vector[N] b_pr;          # go bias 
  vector[N] rho_pr;        # rho, inv temp 
}
transformed parameters{
  vector<lower=0,upper=1>[N] xi;
  vector<lower=0,upper=1>[N] ep;
  vector[N] b; 
  vector<lower=0>[N] rho;
  
  for (i in 1:N) {
    xi[i]  = Phi_approx( mu_p[1] + sigma[1] * xi_pr[i] );
    ep[i]  = Phi_approx( mu_p[2] + sigma[2] * ep_pr[i] );
  }
  b   = mu_p[3] + sigma[3] * b_pr; # vectorization
  rho = exp( mu_p[4] + sigma[4] * rho_pr );
}
model {  
  # gng_m2: RW + noise + bias model in Guitart-Masip et al 2012
  # hyper parameters
  sigma ~ cauchy(0, 5.0);
  
  # individual parameters w/ Matt trick
  xi_pr  ~ normal(0, 1.0);   
  ep_pr  ~ normal(0, 1.0);   
  b_pr   ~ normal(0, 1.0);   
  rho_pr ~ normal(0, 1.0);
  
  # for (i in 1:N) {
  #   vector[K] wv_g;  # action wegith for go
  #   vector[K] wv_ng; # action wegith for nogo
  #   vector[K] qv_g;  # Q value for go
  #   vector[K] qv_ng; # Q value for nogo
  #   vector[K] pGo;   # prob of go (press) 
  #   
  #   wv_g  = initV;
  #   wv_ng = initV;
  #   qv_g  = initV;
  #   qv_ng = initV;
  #   
  #   for (t in 1:Tsubj[i])  {
  #     wv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + b[i];
  #     wv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ];  # qv_ng is always equal to wv_ng (regardless of action)      
  #     pGo[ cue[i,t] ]   = inv_logit( wv_g[ cue[i,t] ] - wv_ng[ cue[i,t] ] ); 
  #     pGo[ cue[i,t] ]   = pGo[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
  #     print(pGo[ cue[i,t] ])
  #     #pressed_r[i,t] ~ bernoulli( pGo[ cue[i,t] ] );
  #     
  #     # update action values
  #     
  #     qv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + ep[i] * (rho[i] * outcome_r[i,t] - qv_g[ cue[i,t] ]);
  #     qv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ] + ep[i] * (rho[i] * outcome_l[i,t] - qv_ng[ cue[i,t] ]);  
  #       
  #   } # end of t loop
  # } # end of i loop
}
generated quantities {
  real<lower=0, upper=1> mu_xi;
  real<lower=0, upper=1> mu_ep;
  real mu_b; 
  real<lower=0> mu_rho;
  real log_lik[N];
  real Qgo[N, T];
  real Qnogo[N, T];
  real Wgo[N, T];
  real Wnogo[N, T];
  int<lower=0, upper=1> pressed_r[N, T];

  mu_xi  = Phi_approx(mu_p[1]);
  mu_ep  = Phi_approx(mu_p[2]);
  mu_b   = mu_p[3];
  mu_rho = exp(mu_p[4]); 
  
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
        wv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + b[i];
        wv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ];  # qv_ng is always equal to wv_ng (regardless of action)      
        pGo[ cue[i,t] ]   = inv_logit( wv_g[ cue[i,t] ] - wv_ng[ cue[i,t] ] ); 
        pGo[ cue[i,t] ]   = pGo[ cue[i,t] ] * (1 - xi[i]) + xi[i]/2;  # noise
        pressed_r[i,t] = bernoulli_rng( pGo[ cue[i,t] ] );
        #log_lik[i] = log_lik[i] + bernoulli_lpmf( pressed_r[i,t] | pGo[ cue[i,t] ] );
        
        # Model regressors --> store values before being updated
        Qgo[i,t]   = qv_g[ cue[i,t] ];
        Qnogo[i,t] = qv_ng[ cue[i,t] ];
        Wgo[i,t]   = wv_g[ cue[i,t] ];
        Wnogo[i,t] = wv_ng[ cue[i,t] ];
        
        # update action values
        if (pressed_r[i,t]) { # update go value 
          qv_g[ cue[i,t] ]  = qv_g[ cue[i,t] ] + ep[i] * (rho[i] * outcome_r[i,t] - qv_g[ cue[i,t] ]);
        } else { # update no-go value  
          qv_ng[ cue[i,t] ] = qv_ng[ cue[i,t] ] + ep[i] * (rho[i] * outcome_l[i,t] - qv_ng[ cue[i,t] ]);  
        }  
      } # end of t loop
    } # end of i loop
  } # end of local section
}'


# data {
#   int<lower=1> N;
#   int<lower=1> T;
#   int<lower=1> K;
#   int<lower=1, upper=T> Tsubj[N];
#   real outcome_l[N, T];
#   real outcome_r[N, T];
#   #int<lower=0, upper=1> pressed_r[N, T];
#   int<lower=1, upper=K> cue[N, T];
# }

p_right <- data.frame(expand.grid(cue=1:6, reward = c(1,5)), pcorrect_if_pressed_r=c(rep(.2,3), rep(.8,3)))

Trials <- p_right[sample(1:dim(p_right)[1], size = 8*48, replace = T),]
Trials$crct_if_right <- rbinom(dim(Trials)[1], size = 1, prob = Trials$pcorrect_if_pressed_r)
Trials$outcome_r <- Trials$crct_if_right*Trials$reward
Trials$outcome_l <- (1-Trials$crct_if_right)*Trials$reward

N <- 2
T <- dim(Trials)[1]
rw_data <- list(N = N,
                T = T,
                K = length(unique(Trials$cue)),
                Tsubj = rep(T, N),
                outcome_l = matrix(rep(Trials$outcome_l, each = N), nrow = N),
                outcome_r = matrix(rep(Trials$outcome_r, each = N), nrow = N),
                cue = matrix(rep(Trials$cue, each =N), nrow = N))

rw_gen_fit <- stan(model_code = rw_model, data = rw_data,  iter = 1000, chains = 8)

somedats <- rstan::extract(rw_gen_fit)

dim(somedats$pressed_r)
dim(somedats$Qgo)

library(tidyverse)

qplot(somedats$Qnogo[1,1,], 1:384)

press_right <- dplyr::bind_rows(lapply(1:4000,function(x){
  tmpdf <- data.frame(cue = Trials$cue, 
                      s1 = somedats$pressed_r[x,1,],
                      s2 = somedats$pressed_r[x,2,])
  tmpdf_sum <- tmpdf %>% 
    gather(subject, pr, -cue) %>%
    group_by(cue, subject) %>%
    summarize(prop_pr = sum(pr)/n())
  tmpdf_sum$iter <- x
  tmpdf_sum
  }))

head(press_right)

lr <- data.frame(s1 = somedats$ep_pr[,1],
                 s2 = somedats$ep_pr[,1],
                 iter = 1:dim(somedats$ep_pr)[1]) %>%
  gather(subject, ep_pr, -iter)

pr_lr <- left_join(press_right, lr)
head(pr_lr)

ggplot(pr_lr, aes(x = ep_pr, y = prop_pr)) +
  geom_point(alpha=.1) +
  geom_smooth(method = 'lm') +
  facet_grid(cue~subject)


press_right_l <- dplyr::bind_rows(lapply(1:4000,function(x){
  tmpdf <- data.frame(cue = Trials$cue, 
                      s1 = somedats$pressed_r[x,1,],
                      s2 = somedats$pressed_r[x,2,],
                      t = 1:length(Trials$cue))
  tmpdf$iter <- x
  tmpdf
})) %>% gather(subject, pr, -cue, -t,-iter)

head(press_right_l)

iters <- sample(unique(press_right_l$iter), size = 100)

ggplot(filter(press_right_l, iter %in% iters), 
       aes(x = t, y = pr, group = iter)) +
  geom_point() +
  geom_line(stat = 'smooth', method = 'loess', alpha = .02) +
  facet_grid(cue~subject)

  