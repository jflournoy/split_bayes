source('load_split_data.r')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

htRW_m1_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg.stan', 
                    data = ht_stan_data, 
                    iter = 10000, chains = 28, control = list(adapt_delta = 0.90), init = genInitListFunc(dl_stan_data$N))
saveRDS(htRW_m1_fit, htRW_m1_fname)