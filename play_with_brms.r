library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('load_split_data.r')
help(brms)

head(stanSplitDF) %>% knitr::kable()



fit <- brm(pressed_r ~ 1 + condition*cue + (1 | id), 
           data = mutate(stanSplitDF, cue = factor(cue)), 
           family = 'bernoulli',
           save_model = 'brm_test.stan',
           save_dso = T,
           iter = 2000,
           chains =4)
summary(fit)
