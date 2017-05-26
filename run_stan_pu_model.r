source('load_split_data.r')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

puRW_m1_fit <- stan(file = file.path(stan_model_dir, 'gng_m1_reg.stan'), 
                    data = pu_stan_data, 
                    iter = 7500, warmup = 5000, chains = 28, control = list(adapt_delta = 0.95), init = genInitListFunc(dl_stan_data$N))
saveRDS(puRW_m1_fit, puRW_m1_fname)
