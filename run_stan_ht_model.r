source('load_split_data.r')

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fname_full <- paste0(htRW_m1_fname, as.character(floor(as.numeric(Sys.time())*1000)), ".RDS")

htRW_m1_fit <- stan(file = file.path(stan_model_dir, 'gng_m1_reg.stan'), 
                    data = ht_stan_data, 
                    iter = 7500, warmup = 5000, chains = 1, control = list(adapt_delta = 0.95), init = genInitListFunc(dl_stan_data$N))
saveRDS(htRW_m1_fit, fname_full)
