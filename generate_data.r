#----Generate RW Data----
#' # Generate data like RW model
#' 
#' 

p_right <- data.frame(expand.grid(cue=1:6, reward = c(1,5)), pcorrect_if_pressed_r=c(rep(.2,3), rep(.8,3)))

Trials <- p_right[sample(1:dim(p_right)[1], size = 8*48, replace = T),]
Trials$crct_if_right <- rbinom(dim(Trials)[1], size = 1, prob = Trials$pcorrect_if_pressed_r)
Trials$outcome_r <- Trials$crct_if_right*Trials$reward
Trials$outcome_l <- (1-Trials$crct_if_right)*Trials$reward

#'
#' ## Play with the model
#' 

library(tidyverse)

inv_logit <- function(x) exp(x)/(1+exp(x))
Phi_approx <- function(x) pnorm(x)

rw_strategy <- function(trialdf, mu_p){
  xi <- Phi_approx( mu_p[1])# + sigma[1] * xi_pr[i] )
  ep <- Phi_approx( mu_p[2])# + sigma[2] * ep_pr[i] )
  b <- mu_p[3]# + sigma[3] * b_pr; # vectorization
  rho <- exp( mu_p[4])# + sigma[4] * rho_pr );
  
  K <- length(unique(trialdf$cue))
  Tsubj <- dim(trialdf)[1]
  wv_g <- c(rep(0, K))  # action wegith for go
  wv_ng <- c(rep(0, K)) # action wegith for nogo
  qv_g <- c(rep(0, K))  # Q value for go
  qv_ng <- c(rep(0, K)) # Q value for nogo
  pGo <- c(rep(0, K))   # prob of go (press)
  
  trialdf$pressed_r <- NA
  trialdf$Qgo       <- NA
  trialdf$Qnogo     <- NA
  trialdf$Wgo       <- NA
  trialdf$Wnogo     <- NA
  trialdf$pGo       <- NA
  trialdf$outcome   <- NA
  
  for (t in 1:Tsubj)  {
    wv_g[ trialdf$cue[t] ] <- qv_g[ trialdf$cue[t] ] + b
    wv_ng[ trialdf$cue[t] ] <-  qv_ng[ trialdf$cue[t] ]  # qv_ng is always equal to wv_ng (regardless of action)
    pGo[ trialdf$cue[t] ]   = inv_logit( wv_g[ trialdf$cue[t] ] - wv_ng[ trialdf$cue[t] ] )
    pGo[ trialdf$cue[t] ]   = pGo[ trialdf$cue[t] ] * (1 - xi) + xi/2;  # noise

    trialdf$pressed_r[t] <- rbinom(n = 1, size = 1, prob = , pGo[ trialdf$cue[t] ]);
    
    trialdf$Qgo[t]   <- qv_g[ trialdf$cue[t] ];
    trialdf$Qnogo[t] <- qv_ng[ trialdf$cue[t] ];
    trialdf$Wgo[t]   <- wv_g[ trialdf$cue[t] ];
    trialdf$Wnogo[t] <- wv_ng[ trialdf$cue[t] ];
    trialdf$pGo[t]   <- pGo[ trialdf$cue[t] ];
    
    # update action values
    if(trialdf$pressed_r[t] == 1){
      qv_g[ trialdf$cue[t] ]  <- qv_g[ trialdf$cue[t] ] + ep * (rho * trialdf$outcome_r[t] - qv_g[ trialdf$cue[t] ]);
      trialdf$outcome[t] <- trialdf$outcome_r[t]
    } else {
      qv_ng[ trialdf$cue[t] ] <- qv_ng[ trialdf$cue[t] ] + ep * (rho * trialdf$outcome_l[t] - qv_ng[ trialdf$cue[t] ]); 
      trialdf$outcome[t] <- trialdf$outcome_l[t]
    }
  } # end of t loop
  return(trialdf)
}

library(gridExtra)
library(cowplot)

plot_RW_run <- function(trials, mu_p){
  single_run <- rw_strategy(trialdf = Trials,
                            mu_p = mu_p)
  
  aplot <- single_run %>%
    filter(cue %in% c(1,4)) %>%
    mutate(cue = factor(cue)) %>%
    group_by(cue) %>%
    mutate(t = 1:n(), last_outcome = as.numeric( ifelse(lag(pressed_r) == 1 & lag(outcome) > 0, 1,
                                                        ifelse(lag(pressed_r) == 1 & lag(outcome) == 0, .1,
                                                               ifelse(lag(pressed_r) == 0 & lag(outcome) > 0, 0,
                                                                      .9)))),
           last_press = lag(pressed_r)) %>%
    ggplot(aes(x = t, y = pGo)) + 
    geom_segment(aes(xend = t, yend = last_outcome), alpha = .1, color = 'black') + 
    geom_point(aes(y = last_outcome, shape = factor(last_press))) + 
    scale_shape_manual(values = c(25,24))+
    scale_y_continuous(breaks = c(0,1), labels = c('left', 'right'))+
    geom_point() + 
    geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, k = 8,  bs = "cs")) + 
    facet_wrap(pcorrect_if_pressed_r~cue, nrow = 2)+
    theme(panel.background = element_blank(), strip.text = element_blank())
  return(aplot)
}
plot_000 <- plot_RW_run(trials = Trials,
                        mu_p = c(xi = 0, ep = 0, b = 0, rho = 0))
plot_000

plot_est <- plot_RW_run(trials = Trials,
                        mu_p = c(xi = qnorm(.1), ep = qnorm(0.045), b = 0, rho = log(1.4)))


plot_higher_ep <- plot_RW_run(trials = Trials,
                              mu_p = c(xi = qnorm(.1), ep = qnorm(0.15), b = 0, rho = log(1.4)))
plot_higher_ep_hr <- plot_RW_run(trials = Trials,
                              mu_p = c(xi = qnorm(.1), ep = qnorm(0.15), b = 0, rho = log(3)))
plot_higher_ep_lr <- plot_RW_run(trials = Trials,
                                 mu_p = c(xi = qnorm(.1), ep = qnorm(0.15), b = 0, rho = log(.5)))

plot_lower_ep <- plot_RW_run(trials = Trials,
                             mu_p = c(xi = qnorm(.1), ep = qnorm(0.03), b = 0, rho = log(1.4)))
plot_lower_ep_hr <- plot_RW_run(trials = Trials,
                             mu_p = c(xi = qnorm(.1), ep = qnorm(0.03), b = 0, rho = log(3)))
plot_lower_ep_lr <- plot_RW_run(trials = Trials,
                             mu_p = c(xi = qnorm(.1), ep = qnorm(0.03), b = 0, rho = log(.5)))


plot_higher_rho <- plot_RW_run(trials = Trials,
                               mu_p = c(xi = qnorm(.1), ep = qnorm(0.045), b = 0, rho = log(3)))
plot_lower_rho <- plot_RW_run(trials = Trials,
                              mu_p = c(xi = qnorm(.1), ep = qnorm(0.045), b = 0, rho = log(.5)))

plot_grid(plot_higher_ep_lr, plot_higher_ep, plot_higher_ep_hr,
          plot_lower_rho, plot_est, plot_higher_rho,
          plot_lower_ep_lr, plot_lower_ep, plot_lower_ep_hr,
          nrow = 3, labels = c('he_lr', 'he', 'he_hr', 'lr', 'est', 'hr', 'le_lr', 'le', 'le_hr'))

#----Stan Generate RW Data----
#'
#' ## Stan generate using these params
#' 

rw_model <- '
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> K;
  int<lower=1, upper=T> Tsubj[N];
  real outcome_r[N, T];
  real outcome_l[N, T];
  #int<lower=0, upper=1> pressed_r[N, T];
  int<lower=1, upper=K> cue[N, T];
  vector[4] mu_p; 
}
transformed data {
  vector[K] initV;
  initV  = rep_vector(0.0, K);
}
parameters {
  # declare as vectors for vectorizing
  # vector[4] mu_p;  
  #vector<lower=0>[4] sigma; 
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
    xi[i]  = Phi_approx( mu_p[1] + xi_pr[i] );
    ep[i]  = Phi_approx( mu_p[2] + ep_pr[i] );
  }
  b   = mu_p[3] +  b_pr; # vectorization
  rho = exp( mu_p[4] +  rho_pr );
}
model {  
  # gng_m2: RW + noise + bias model in Guitart-Masip et al 2012
  # hyper parameters
  
  # individual parameters w/ Matt trick
  xi_pr  ~ normal(0, 0.25);   
  ep_pr  ~ normal(0, 0.25);   
  b_pr   ~ normal(0, 0.25);   
  rho_pr ~ normal(0, 0.25);
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
  real Pgo[N, T];
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
        Pgo[i,t]   = pGo[ cue[i,t] ];
        
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

N <- 2
T <- dim(Trials)[1]
rw_data <- list(N = N,
                T = T,
                K = length(unique(Trials$cue)),
                Tsubj = rep(T, N),
                outcome_l = matrix(rep(Trials$outcome_l, each = N), nrow = N),
                outcome_r = matrix(rep(Trials$outcome_r, each = N), nrow = N),
                cue = matrix(rep(Trials$cue, each =N), nrow = N),
                mu_p = c(-1, -2, 0, .4))

rw_gen_fname <- file.path('/data/jflournoy/split/bayes/', 'rw_gen_stan.RDS')
if(file.exists(rw_gen_fname)){
  rw_gen_fit <- readRDS(rw_gen_fname)
} else {
  rw_gen_fit <- stan(model_code = rw_model, data = rw_data,  iter = 1000, chains = 8)
  saveRDS(rw_gen_fit, rw_gen_fname)
}

rw_gen_fit@model_pars

stan_trace(rw_gen_fit, pars = c('xi', 'ep', 'b','rho'))
stan_dens(rw_gen_fit, pars = c('xi', 'ep', 'b','rho'))

stan_dens(rw_gen_fit, pars = c('mu_rho', 'rho'))

somedats <- rstan::extract(rw_gen_fit)

library(tidyverse)

pgo_l <- dplyr::bind_rows(lapply(1:4000,function(x){
  tmpdf <- data.frame(cue = Trials$cue, 
                      s1 = somedats$Pgo[x,1,],
                      s2 = somedats$Pgo[x,2,],
                      t = 1:length(Trials$cue))
  tmpdf$iter <- x
  tmpdf
})) %>% gather(subject, pr, -cue, -t,-iter)

#+cache=T
pgo_l %>%
  filter(iter %in% base::sample(unique(iter), size = 500)) %>%
  ggplot(aes(x = t, y = pr, group = iter)) + 
  geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, k = 8,  bs = "cs"), alpha=.25) + 
  facet_wrap(subject~cue, nrow = 2)+
  theme(panel.background = element_blank())

#----test parameter recovert----
#'
#' # Can we recover parameters
#' 

Trials_ct <- Trials %>%
  mutate(cuetype = c(1:3, 1:3)[cue])

head(Trials_ct)

pressed_right <- dplyr::bind_rows(lapply(1:4000,function(x){
  snum <- sample(c(1,2), size = 1)
  pr <- somedats$pressed_r[x,snum,]
  outcome <- ifelse(somedats$pressed_r[x,snum,], Trials_ct$outcome_r, Trials_ct$outcome_l)
  tmpdf <- data.frame(cue = Trials$cue,
                      pr = pr,
                      outcome = outcome,
                      t = 1:length(Trials$cue))
  tmpdf$Subj <- x
  tmpdf
}))

head(pressed_right)
dim(pressed_right)

subject_nums <- sample(1:4000, size = 100)

outcome <- pressed_right %>%
  select(t, outcome, Subj) %>%
  spread(t, outcome) %>%
  slice(subject_nums) %>%
  select(-Subj) %>%
  as.matrix
head(outcome)
dim(outcome)

pressed_r <- pressed_right %>%
  select(t, pr, Subj) %>%
  spread(t, pr) %>%
  slice(subject_nums) %>%
  select(-Subj) %>%
  as.matrix
head(pressed_r)
dim(outcome)

#----compare timing with different "vectorization"----
#'
#' ## First, compare timing for vectorized versus unvectorized hyperparameters
#' 

N <- dim(outcome)[1]
T <- dim(outcome)[2]
rw_hb_test_data <- list(N = N,
                     T = T,
                     K = length(unique(Trials_ct$cue)),
                     Tsubj = rep(T, N),
                     outcome = outcome,
                     pressed = pressed_r,
                     cue = matrix(rep(Trials_ct$cue, each = N), nrow = N))

rw_hb_test_fname <- file.path('/data/jflournoy/split/bayes/', 'rw_hb_test_stan.RDS')
if(file.exists(rw_hb_test_fname)){
  rw_hb_test_fit <- readRDS(rw_hb_test_fname)
} else {
  rw_hb_test_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg.stan', data = rw_hb_test_data,  iter = 10, chains = 1)
  saveRDS(rw_hb_test_fit, rw_hb_test_fname)
}

print(get_elapsed_time(rw_hb_test_fit))

rw_hb_uv_test_fname <- file.path('/data/jflournoy/split/bayes/', 'rw_hb_uv_test_stan.RDS')
if(file.exists(rw_hb_uv_test_fname)){
  rw_hb_uv_test_fit <- readRDS(rw_hb_uv_test_fname)
} else {
  rw_hb_uv_test_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg_unvecparam.stan', data = rw_hb_test_data,  iter = 10, chains = 1)
  saveRDS(rw_hb_uv_test_fit, rw_hb_uv_test_fname)
}

print(get_elapsed_time(rw_hb_uv_test_fit))

#----recover parameter estimates from hBayesDM model----
#'
#' ## Get model estiamtes from the basic RW model
#' 

rw_hb_uv_test_1000_fname <- file.path('/data/jflournoy/split/bayes/', 'rw_hb_uv_test_1000_stan.RDS')
if(file.exists(rw_hb_uv_test_1000_fname)){
  rw_hb_uv_test_1000_fit <- readRDS(rw_hb_uv_test_1000_fname)
} else {
  rw_hb_uv_test_1000_fit <- stan(fit = rw_hb_uv_test_fit, data = rw_hb_test_data,  iter = 1000, chains = 1)
  saveRDS(rw_hb_uv_test_1000_fit, rw_hb_uv_test_1000_fname)
}

stan_dens(rw_hb_uv_test_1000_fit, pars = c('mu_xi', 'mu_ep', 'mu_rho'))

#+rho plot, fig.width= 20, fig.height = 20
stan_dens(rw_hb_uv_test_1000_fit, pars = c('ep'))

#----Use very simple model to check cholesky_corr reparameterization----
#' 
#' ## Check reparameterization in simple model
#'
#' We want to add a mean for each cue type, but it wasn't playing nice in the big model, so let's just make sure we're 
#' setting it up correctly.
#'   

# data {
#   int<lower=1> N; #Number of subjects
#   int<lower=1> T; #Max number of trials
#   int<lower=1> K; #Number of distinct cues
#   int<lower=1> D; #Number of groups of cues
#   int<lower=1, upper=T> Tsubj[N]; #number of trials per subject
#   int<lower=0, upper=1> pressed_r[N, T]; #Subject x Trial matrix of button presses (L = 0, R = 1)
#   int<lower=1, upper=K> cue[N, T]; #Subject x Trial matrix of which cue was presented
#   int<lower=1, upper=D> cuetype[N, T]; #Subject x Trial matrix of which type of cue was presented
# }

library(MASS)
mu_cuetype <- c(-1,0,1)
Sigma_cuetype <- matrix(c(1,.1,.2,  .1,1,-.1, .2,-.1,1), ncol = 3)

N <- 100
T <- 150
K <- 6
D <- 3
Tsubj <- rep(T, N)
cues <- data.frame(cue = 1:6, cuetype = rep(1:3, 2))
trial_structure <- cues[sample(1:6, size = T, replace = T), ]
cue <- matrix(rep(trial_structure$cue, each = N), nrow = N)
cuetype <- matrix(rep(trial_structure$cuetype, each = N), nrow = N)

mu_N <-bind_rows(lapply(as.list(1:N), function(x) as.data.frame(t(mvrnorm(n = 1, mu = mu_cuetype, Sigma = Sigma_cuetype))) ))

bin_outcome <- matrix(nrow = N, ncol = T)

for(n in 1:N){
  for(t in 1:T){
    ct <- cuetype[n,t]
    pR <- inv_logit(Phi_approx(mu_N[n,ct]))
    bin_outcome[n,t] <- rbinom(n = 1, size = 1, prob = pR)
  }
}

logit <- function(x) log(x/(1-x))

qnorm(logit(mean(bin_outcome[,trial_structure$cuetype==1])))
qnorm(logit(mean(bin_outcome[,trial_structure$cuetype==2])))
qnorm(logit(mean(bin_outcome[,trial_structure$cuetype==3])))

test_reparam_data <- list(N = N, T = T, K = K, D = D, Tsubj = Tsubj,  pressed_r = bin_outcome, cue = cue, cuetype = cuetype)

test_reparam_fname <- file.path('/data/jflournoy/split/bayes/', 'test_reparam_chol_stan.RDS')
if(file.exists(test_reparam_fname)){
  test_reparam_fit <- readRDS(test_reparam_fname)
} else {
  test_reparam_fit <- stan(file = '~/code_new/split_bayes/test_reparam_chol.stan', data = test_reparam_data,  iter = 1000, chains = 6)
  saveRDS(test_reparam_fit, test_reparam_fname)
}

test_reparam_fit@model_pars
test_reparam_fit@par_dims

stan_ac(test_reparam_fit, pars = c('mu_ep', 'L_Omega_ep', 'tau_ep'))
stan_trace(test_reparam_fit, pars = c('mu_ep', 'L_Omega_ep', 'tau_ep'))
stan_dens(test_reparam_fit, pars = c('mu_ep','Sigma_ep'))
stan_plot(test_reparam_fit, pars = c('mu_ep'))
stan_plot(test_reparam_fit, pars = c('L_Omega_ep'))

stan_plot(test_reparam_fit, pars = c('some_sig'))
stan_plot(test_reparam_fit, pars = c('indiv_var'))
stan_plot(test_reparam_fit, pars = c('ep_raw'))
stan_plot(test_reparam_fit, pars = c('other_sig'))

eps <- rstan::extract(test_reparam_fit, pars = 'ep')

qnorm(mean(eps$ep[,,1]))
qnorm(mean(eps$ep[,,2]))
qnorm(mean(eps$ep[,,3]))

#----Finally test on modified gng_m1----
#'
#' ## Finally test on the modified gng_m1
#' 

outcome <- pressed_right %>%
  dplyr::select(t, outcome, Subj) %>%
  spread(t, outcome) %>%
  slice(subject_nums) %>%
  dplyr::select(-Subj) %>%
  as.matrix
head(outcome)
dim(outcome)


pressed_r <- pressed_right %>%
  dplyr::select(t, pr, Subj) %>%
  spread(t, pr) %>%
  slice(subject_nums) %>%
  dplyr::select(-Subj) %>%
  as.matrix
head(pressed_r)
dim(outcome)

# data {
#   int<lower=1> N; #Number of subjects
#   int<lower=1> T; #Max number of trials
#   int<lower=1> K; #Number of distinct cues
#   int<lower=1> D; #Number of groups of cues
#   int<lower=1, upper=T> Tsubj[N]; #number of trials per subject
#   real outcome[N, T]; #Subject x Trial matrix of outcomes
#   int<lower=0, upper=1> pressed_r[N, T]; #Subject x Trial matrix of button presses (L = 0, R = 1)
#   int<lower=1, upper=K> cue[N, T]; #Subject x Trial matrix of which cue was presented
#   int<lower=1, upper=D> cuetype[N, T]; #Subject x Trial matrix of which type of cue was presented
# }

N <- dim(outcome)[1]
T <- dim(outcome)[2]
rw_test_data <- list(N = N,
                     T = T,
                     K = length(unique(Trials_ct$cue)),
                     D = length(unique(Trials_ct$cuetype)),
                     Tsubj = rep(T, N),
                     outcome = outcome,
                     pressed_r = pressed_r,
                     cue = matrix(rep(Trials_ct$cue, each = N), nrow = N),
                     cuetype = matrix(rep(Trials_ct$cuetype, each = N), nrow = N))

rw_test_fname <- file.path('/data/jflournoy/split/bayes/', 'rw_test_stan.RDS')
if(file.exists(rw_test_fname)){
  rw_test_fit <- readRDS(rw_test_fname)
} else {
  rw_test_fit <- stan(file = '~/code_new/split_bayes/split_m1_repar_reg.stan', 
                      data = rw_test_data,  iter = 10, chains = 1)
  saveRDS(rw_test_fit, rw_test_fname)
}

pairs(rw_test_fit, pars = c('mu_ep'))

#'   Warning messages:
#'   1: There were 3000 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
#'   http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
#'   2: Examine the pairs() plot to diagnose sampling problems

rw_test_fit@model_pars
rw_test_fit@par_dims
stan_trace(rw_test_fit, pars = c('mu_ep'))

pairs(rw_test_fit, pars = c('mu_ep', 'lp__'))
pairs(rw_test_fit, pars = c('mu_ep', 'lp__'))
