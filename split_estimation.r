#' ---
#' Title: Estimates from TDS-II SPLiT task
#' Author: John C Flournoy
#' ---
#' 
#+set up data extraction,echo=F,warning=F,error=F----

library(tidyverse)
library(stringr)
library(knitr)
library(rstan)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

age_data <- read_csv('/data/jflournoy/split/age_gender_iq.csv')

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

excl_ids <- c(110) # 110 didn't play the same game
root_dir <- '/data/jflournoy/split/' # This would be '/Volumes/' if you've mounted locally
raw_data_dir <- paste0(root_dir, 'experiment_data/')
processed_data_dir <-  paste0(root_dir, 'processed/')
#previous iterations of this task included 'BannasOranges' and we don't want those trials
conditions <- c('HungryThirsty', 'DatingLooking', 'PopularUnpopular')
#We will use this later to separate the metacog questions from trial behavior
metacog_condition <- 'Confidence'
file_match_pattern <- 'split.*csv'

#This turns estimates on the logit scale back into probabilities
invlogit <- function(x) {
  exp(x)/(exp(x)+1)
}

files <- data_frame(file=dir(path=raw_data_dir, pattern=file_match_pattern)) 

col_spec <- cols(
  rt = col_integer(),
  key_press = col_integer(),
  trial_type = col_character(),
  trial_index = col_integer(),
  time_elapsed = col_integer(),
  internal_node_id = col_character(),
  subject = col_character(),
  responses = col_character(),
  trial_id = col_character(),
  text = col_character(),
  block_duration = col_integer(),
  timing_post_trial = col_integer(),
  view_history = col_character(),
  correct = col_character(),
  stimulus = col_character(),
  correct_response = col_integer(),
  possible_responses = col_number(),
  stim_duration = col_integer(),
  feedback_duration = col_integer(),
  exp_stage = col_character(),
  image = col_character(),
  context = col_character(),
  condition = col_character(),
  optimal_response = col_integer(),
  stim_gender = col_character(),
  stimindex = col_integer(),
  feedback = col_character(),
  reward_possible = col_integer(),
  stim_chosen = col_integer(),
  exp_id = col_character(),
  qnum = col_integer(),
  page_num = col_integer(),
  trial_num = col_integer(),
  stim_question = col_character(),
  stim_response = col_character(),
  score_response = col_integer(),
  response_range = col_number(),
  times_viewed = col_integer(),
  credit_var = col_character(),
  total_points = col_character()
)

splitDF_raw <- files %>% 
  tidyr::extract(file, c('filename', 'id', 'timestamp'), 
          '(split-([0-9]{3}|bad_pid)[-_]([0-9]+)\\.csv)') %>% 
  filter(id < 500 & !id %in% excl_ids) %>%
  group_by(filename, id, timestamp) %>%
  do({
    # print(.$id[[1]])
    dat=read_csv(paste0(raw_data_dir, .$filename[[1]]), col_types = col_spec)
    # data_frame(dat = list(dat))
    dat
    })

splitDF_ <- splitDF_raw %>%
  filter((trial_id == 'stim' & (key_press>0 | key_press==-1)) | (trial_id == 'metacog')) %>%
  tidyr::extract(image, c('stim_image', 'sex'), '.*/(([fma])_f42887_e_[0-9]{3}\\.png)') %>%
  mutate(correcttrue=ifelse(key_press==-1, NA, as.numeric(correct=='true')),
         outcome=ifelse(feedback, reward_possible, 0),
         pressed_r=c(`37`=0,`39`=1)[as.character(key_press)]) %>%
  rowwise %>%
  mutate(condition=paste(sort(unlist(str_split(context, '_'))), collapse=''), 
         condition=ifelse(trial_id=='metacog', 'Confidence', condition),
         correcttrue=ifelse(trial_id=='metacog', 
                            as.numeric(score_response), correcttrue)) %>%
  filter(condition %in% c(conditions, metacog_condition)) %>%
  group_by(id, condition) %>%
  mutate(condition_trial_index=1:n(),
         block=((condition_trial_index-1) %/% 16) +1) %>%
  ungroup() %>%
  select(id, correcttrue, pressed_r, outcome, condition, condition_trial_index, sex, stim_image, trial_index, block)


metacogDF <- splitDF_ %>% 
  select(-sex, -stim_image) %>%
  filter(condition == metacog_condition,
         !is.na(correcttrue)) %>%
  mutate(block=(condition_trial_index+1)/2) %>% 
  spread(condition, correcttrue) %>%
  rename(confidence=Confidence) %>%
  select(-condition_trial_index, -trial_index)

splitDF<- splitDF_ %>% filter(condition != metacog_condition) %>%
  left_join(metacogDF) %>%
  mutate(scaled_index=scale(condition_trial_index),
         condition=factor(condition, levels=conditions, labels=abbreviate(conditions))) %>%
  arrange(id, trial_index)

splitDF %>% filter(!is.na(correcttrue)) %>%
  group_by(id) %>% 
  distinct(condition) %>%
  summarise(n_conds=n()) %>%
  kable

splitDF %>%  filter(!is.na(correcttrue)) %>%
  group_by(id, condition, block) %>%
  summarise(n_trials=n(),
            pcor=mean(correcttrue, na.rm=T),
            confidence=unique(confidence)) %>%
  #   select(-n_trials) %>%
  gather(key, value, n_trials, pcor, confidence) %>%
  unite(cond_stat, condition, key) %>%
  spread(cond_stat, value) %>%
  kable(digits=2)

#'
#' # Bayesin estimation of Rescorla-Wagner process
#' 
#' At the moment, there is a separate model for each For now, running in separate models per condition
#'
#+make data nice for stan,echo=T----

stanSplitDF <- splitDF %>% 
  filter(!is.na(pressed_r)) %>%
  group_by(id, condition) %>%
  mutate(condition_index = 1:n(),
         cue = as.numeric(factor(sex, levels = c('m','f'), labels = c('m','f'))))
  
condition_sum <- stanSplitDF  %>%
  group_by(id, condition) %>%
  summarize(K=length(unique(stim_image)),
            Tsubj=n()) %>%
  group_by(condition) %>%
  mutate(N = length(unique(id)),
         T = max(Tsubj))

outcomeDFs <- stanSplitDF %>%
  select(id,outcome,condition_index) %>%
  group_by(condition) %>%
  do({
    data_frame( condDF = list( select(spread(ungroup(.), condition_index, outcome, fill = 1), -id, -condition) ) )
  })

pressedrDFs <- stanSplitDF %>%
  select(id,pressed_r,condition_index) %>%
  group_by(condition) %>%
  do({
    data_frame( condDF = list( select(spread(ungroup(.), condition_index, pressed_r, fill = 1), -id, -condition) ) )
  })

cueDFs <- stanSplitDF %>%
  select(id,cue,condition_index) %>%
  group_by(condition) %>%
  do({
    data_frame( condDF = list( select(spread(ungroup(.), condition_index, cue, fill = 1), -id, -condition) ) )
  })

stanAges <- left_join(distinct(ungroup(stanSplitDF),id), mutate(age_data,id=as.character(`subject-name`))) %>%
  select(id,gender,age,iq)

# data {
#   int<lower=1> N;
#   int<lower=1> T;
#   int<lower=2> K;
#   int<lower=1, upper=T> Tsubj[N];
#   real outcome[N, T];
#   int<lower=0, upper=1> pressed[N, T];
#   int<lower=1, upper=K> cue[N, T];
# }

ht_stan_data <- list(N = unique(condition_sum$N[condition_sum$condition == 'HngT']),
                     T = unique(condition_sum$T[condition_sum$condition == 'HngT']),
                     K = unique(condition_sum$K[condition_sum$condition == 'HngT']),
                     Tsubj = condition_sum$Tsubj[condition_sum$condition == 'HngT'],
                     outcome = as.matrix(outcomeDFs$condDF[outcomeDFs$condition == 'HngT'][[1]]),
                     pressed = as.matrix(pressedrDFs$condDF[outcomeDFs$condition == 'HngT'][[1]]),
                     cue = as.matrix(cueDFs$condDF[outcomeDFs$condition == 'HngT'][[1]]))

dl_stan_data <- list(N = unique(condition_sum$N[condition_sum$condition == 'DtnL']),
                     T = unique(condition_sum$T[condition_sum$condition == 'DtnL']),
                     K = unique(condition_sum$K[condition_sum$condition == 'DtnL']),
                     Tsubj = condition_sum$Tsubj[condition_sum$condition == 'DtnL'],
                     outcome = as.matrix(outcomeDFs$condDF[outcomeDFs$condition == 'DtnL'][[1]]),
                     pressed = as.matrix(pressedrDFs$condDF[outcomeDFs$condition == 'DtnL'][[1]]),
                     cue = as.matrix(cueDFs$condDF[outcomeDFs$condition == 'DtnL'][[1]]))

pu_stan_data <- list(N = unique(condition_sum$N[condition_sum$condition == 'PplU']),
                     T = unique(condition_sum$T[condition_sum$condition == 'PplU']),
                     K = unique(condition_sum$K[condition_sum$condition == 'PplU']),
                     Tsubj = condition_sum$Tsubj[condition_sum$condition == 'PplU'],
                     outcome = as.matrix(outcomeDFs$condDF[outcomeDFs$condition == 'PplU'][[1]]),
                     pressed = as.matrix(pressedrDFs$condDF[outcomeDFs$condition == 'PplU'][[1]]),
                     cue = as.matrix(cueDFs$condDF[outcomeDFs$condition == 'PplU'][[1]]))

genInitListFunc <- function(numSubjs) {
  function(){
    inits_fixed <- c(0.10, 0.20, exp(2.0))
    list(
      mu_p   = c(qnorm(inits_fixed[1]), qnorm(inits_fixed[2]), log(inits_fixed[3])),
      sigma  = c(1.0, 1.0, 1.0),
      xi_pr  = rep(qnorm(inits_fixed[1]), numSubjs),
      ep_pr  = rep(qnorm(inits_fixed[2]), numSubjs),
      rho_pr = rep(log(inits_fixed[3]), numSubjs)
    )
  }
}

htRW_m1_fname <- file.path('/data/jflournoy/split/bayes/', 'htRW_m1_stan.RDS')
if(file.exists(htRW_m1_fname)){
  htRW_m1_fit <- readRDS(htRW_m1_fname)
} else {
  htRW_m1_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg.stan', 
                      data = ht_stan_data, 
                      iter = 1000, chains = 6, control = list(adapt_delta = 0.90), init = genInitListFunc(ht_stan_data$N))
  saveRDS(htRW_m1_fit, htRW_m1_fname)
}

dlRW_m1_fname <- file.path('/data/jflournoy/split/bayes/', 'dlRW_m1_stan.RDS')
if(file.exists(dlRW_m1_fname)){
  dlRW_m1_fit <- readRDS(dlRW_m1_fname)
} else {
  dlRW_m1_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg.stan', 
                      data = dl_stan_data, 
                      iter = 1000, chains = 6, control = list(adapt_delta = 0.90), init = genInitListFunc(dl_stan_data$N))
  saveRDS(dlRW_m1_fit, dlRW_m1_fname)
}

puRW_m1_fname <- file.path('/data/jflournoy/split/bayes/', 'puRW_m1_stan.RDS')
if(file.exists(puRW_m1_fname)){
  puRW_m1_fit <- readRDS(puRW_m1_fname)
} else {
  puRW_m1_fit <- stan(file = '~/code_new/split_bayes/gng_m1_reg.stan', 
                      data = pu_stan_data, 
                      iter = 1000, chains = 6, control = list(adapt_delta = 0.90), init = genInitListFunc(pu_stan_data$N))
  saveRDS(puRW_m1_fit, puRW_m1_fname)
}

stan_trace(htRW_m1_fit, pars = c('mu_xi', 'mu_ep', 'mu_rho'))
stan_trace(htRW_m1_fit, pars = c('mu_p','sigma'))
stan_plot(htRW_m1_fit, pars = c('mu_p'))
stan_plot(htRW_m1_fit, pars = c('mu_xi', 'mu_ep', 'mu_rho'))
stan_plot(htRW_m1_fit, pars = c('mu_ep'))
stan_plot(htRW_m1_fit, pars = c('sigma'))
stan_plot(htRW_m1_fit, pars = c('ep_pr'))
stan_plot(htRW_m1_fit, pars = c('rho_pr'))
stan_plot(htRW_m1_fit, pars = c('xi_pr'))



multiplot(stan_plot(htRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Hungry/Thirsty'),
          stan_plot(dlRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Dating/Looking'),
          stan_plot(puRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Popular/Unpopular'),
          cols=1)
          

multiplot(stan_plot(htRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Hungry/Thirsty')+lims(x=c(.5,3)),
          stan_plot(dlRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Dating/Looking')+lims(x=c(.5,3)),
          stan_plot(puRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Popular/Unpopular')+lims(x=c(.5,3)),
          cols=1)

multiplot(stan_plot(htRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Hungry/Thirsty')+lims(x=c(0,.4)),
          stan_plot(dlRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Dating/Looking')+lims(x=c(0,.4)),
          stan_plot(puRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Popular/Unpopular')+lims(x=c(0,.4)),
          cols=1)

library(bayesplot)

htRW_m1_fit@model_pars

htRW_m1_fit_extract <- extract(htRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)
dlRW_m1_fit_extract <- extract(dlRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)
puRW_m1_fit_extract <- extract(puRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)

samples <- dim(htRW_m1_fit_extract$ep)[1]

hist(htRW_m1_fit_extract$ep)
hist(qnorm(htRW_m1_fit_extract$ep))

hist(htRW_m1_fit_extract$rho)
hist(log(htRW_m1_fit_extract$rho))

cor_ep_rho <- unlist(lapply(1:samples, function(x){
  cor(qnorm(htRW_m1_fit_extract$ep[x,]),log(htRW_m1_fit_extract$rho[x,]))
}))

data_frame(cor=cor_ep_rho) %>%
  do(data.frame(t(quantile(.$cor,probs = c(.025,.1,.5,.9,.975))))) %>%
  ggplot(aes(x = X50., y=1))+
  geom_density(aes(x = cor, y = ..density..), data = data_frame(cor = cor_ep_rho))+
  geom_errorbarh(aes(xmin=X2.5.,xmax=X97.5.), size=1, height=0,col='red')+
  geom_errorbarh(aes(xmin=X10.,xmax=X90.), size=4, height=0,col='red')+
  geom_point(size=5,col='black')+
  theme(panel.background = element_blank())+
  labs(x = 'Correlation of learning rate with inverse temperature', y='')

cor_dist_plot <- function(stan_samples, tocor, labslist = list(title='Correlation Density', x='', y='')){
  nsamples <- dim(stan_samples)[1]
  manycors <- unlist(lapply(1:nsamples, function(x){
    cor(stan_samples[x,],tocor)
  }))
  
  corplot <- data_frame(cor=manycors) %>%
    do(data.frame(t(quantile(.$cor,probs = c(.025,.1,.5,.9,.975))))) %>%
    ggplot(aes(x = X50., y=1))+
    geom_density(aes(x = cor, y = ..density..), data = data_frame(cor = manycors))+
    geom_errorbarh(aes(xmin=X2.5.,xmax=X97.5.), size=1, height=0,col='red')+
    geom_errorbarh(aes(xmin=X10.,xmax=X90.), size=4, height=0,col='red')+
    geom_point(size=5,col='black')+
    theme(panel.background = element_blank())+
    labs(labslist)
  corplot
}

#+Age Plots with ep,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(qnorm(htRW_m1_fit_extract$ep), stanAges$age, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of ep with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(qnorm(dlRW_m1_fit_extract$ep), stanAges$age, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of ep with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(qnorm(puRW_m1_fit_extract$ep), stanAges$age, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of ep with age', y=''))+
            lims(x=c(-.5,.5)))

#+Age Plots with rho,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(log(htRW_m1_fit_extract$rho), stanAges$age, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of rho with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(log(dlRW_m1_fit_extract$rho), stanAges$age, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of rho with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(log(puRW_m1_fit_extract$rho), stanAges$age, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of rho with age', y=''))+
            lims(x=c(-.5,.5)))

#+Age Plots with xi,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$age, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of xi with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$age, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of xi with age', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$age, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of xi with age', y=''))+
            lims(x=c(-.5,.5)))

#+IQ Plots with ep,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(qnorm(htRW_m1_fit_extract$ep), stanAges$iq, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of ep with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(qnorm(dlRW_m1_fit_extract$ep), stanAges$iq, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of ep with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(qnorm(puRW_m1_fit_extract$ep), stanAges$iq, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of ep with iq', y=''))+
            lims(x=c(-.5,.5)))

#+IQ Plots with rho,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(log(htRW_m1_fit_extract$rho), stanAges$iq, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of rho with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(log(dlRW_m1_fit_extract$rho), stanAges$iq, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of rho with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(log(puRW_m1_fit_extract$rho), stanAges$iq, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of rho with iq', y=''))+
            lims(x=c(-.5,.5)))

#+IQ Plots with xi,fig.width=10,fig.height=6----
multiplot(cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$iq, 
                        labslist = list(title = 'Hungry/Thirsty', x = 'Correlation of xi with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$iq, 
                        labslist = list(title = 'Dating/Looking', x = 'Correlation of xi with iq', y=''))+
            lims(x=c(-.5,.5)),
          cor_dist_plot(htRW_m1_fit_extract$xi_pr, stanAges$iq, 
                        labslist = list(title = 'Popular/Unpopular', x = 'Correlation of xi with iq', y=''))+
            lims(x=c(-.5,.5)))

#+Age regressions----
library(wesanderson)
summary(stanAges$age)

regression_predictions <- function(stan_samples, regressor, polyorder = 2){
  predict_regressor <- data_frame(regressor = seq(min(regressor), max(regressor), length.out = 2*length(regressor)))
  nsamples <- dim(htRW_m1_fit_extract$ep)[1]
  manyregs <- bind_rows(lapply(1:nsamples, function(x){
    amod <- lm(stan_samples[x,] ~ 1 + poly(regressor,polyorder))
    preds <- as_data_frame(predict(amod, newdata = predict_regressor, interval = 'conf'))
    preds$iter <- x
    preds
  }))
  manyregs$regressor <- rep(predict_regressor$regressor, nsamples)  
  manyregs
}

reg_pred_plot <- function(regsamples, labslist = list(title='Regression plot', x='', y=''), points_frame = NA, palette = c('red', 'blue', 'grey', 'black', 'black')){
  regsamples_interval <- regsamples %>%
    group_by(regressor) %>%
    summarize(fit50 = median(fit), fit025 = quantile(fit, probs = c(.025)), fit975 = quantile(fit, probs = .974),
              upr975 = quantile(upr, probs = c(.975)), lwr025 = quantile(lwr, probs = .025))
  
  aplot <- regsamples %>%
    ggplot(aes(x = regressor))+
    geom_line(aes(group = iter, y = fit), alpha = .05, color = palette[4])+
    geom_ribbon(aes(ymin = lwr025, ymax = upr975, fill = 'ci95'), alpha = .25, data = regsamples_interval)+
    geom_ribbon(aes(ymin = fit025, ymax = fit975, fill = 'fit95'), alpha = .25, data = regsamples_interval)+
    geom_line(aes(y = fit50, color = 'fit50'), size = 2, data = regsamples_interval)+
    scale_fill_manual(name = 'Intervals', values = c(ci95=palette[3], fit95=palette[2]),
                      labels = c(ci95='95% regression CIs', fit95='95% regression fit'), guide = "legend")+
    scale_color_manual(name = '', values = c(fit50=palette[1]), labels = c(fit50='Median regression fit'), guide = "legend")+
    labs(labslist)+
    theme(panel.background = element_blank())
  if("data.frame" %in% class(points_frame)){
    aplot <- aplot +
      geom_point(aes(y = points), data = points_frame, color = palette[5])
  }
  aplot
}

regpalette <- c(wes_palette('Moonrise1')[1:4], 'black')

#+Age HT regressions,fig.width=10,fig.height=6----
ep_points <- data_frame(points = apply(htRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2 <- regression_predictions(htRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2, points_frame = ep_points, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))

rho_points <- data_frame(points = apply(htRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2 <- regression_predictions(htRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2, points_frame = rho_points, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))

xi_points <- data_frame(points = apply(htRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2 <- regression_predictions(htRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2, points_frame = xi_points, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))

