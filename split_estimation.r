#' ---
#' title: Estimates from TDS-II SPLiT task
#' author: John C Flournoy
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     self_contained: false
#' ---
#' 
#+set up,echo=F,warning=F,error=F,message=F,results='hide'----
library(tidyverse)
library(stringr)
library(knitr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

knitr::read_chunk('load_split_data.r')
knitr::read_chunk('run_stan_ht_model.r')
knitr::read_chunk('run_stan_dl_model.r')
knitr::read_chunk('run_stan_pu_model.r')

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

#Helpers for manipulating 0-1 ranges
logit <- function(x){
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x)/(exp(x)+1)
}

#'
#'
#+Load data,echo=F,warning=F,error=F,message=F,results='hide'----

#'
#'
#'
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
#' # Bayesian estimation of a Rescorla-Wagner learning model
#' 
#' The RW parameters below are estimated separately for each condition. Regressions and correlations are estimated using non-bayesian 
#' estimation methods once per MCMC sample. Estimation of these additional parameters will be brought together under
#' a single model over the coming weeks.
#'
#' The RW model for this task takes advantage of the go-nogo model implemented in the [`hBayesDM`](https://rpubs.com/CCSL/hBayesDM) package -- 
#' specifically the first parameterization, which follows Guitart-Masip et al. (2012).
#' 
#' In the condition-specific models, an action weight is calculated for each of 2 cues on each trial. This action weight governs the probability of 
#' choosing a particular response (i.e., the right arrow key rather than the left) for a particular stimulus. The action weight is modified by
#' a noise parameter $\xi$ that allows for a certain amount of drift back toward chance-responding, and bounds the maximum probability of choosing
#' a particular response. The full model of probability, using the tried-and-true inverse logit is:
#' 
#' $$
#' p(a_t|s_t) = \frac{e^{W(a_t|s_t)}}{1+e^{W(a_t|s_t)}} (1-\xi) + \frac{\xi}{2}
#' $$
#' 
#' where $W$ is the weight for action $a_t = \text{press right-arrow}$ for stimulus $s_t$.
#' 
#' In this particular parameterization, the action weight $W(a,s)$ is just $Q(a,s)$ which is computed for each trial $t$ as follows:
#' 
#' 
#' $$
#' Q_{t}(a_{t},s_{t}) = Q_{t-1}(a_{t},s_{t}) + \epsilon(\rho r_{t} - Q_{t-1}(a_{t},s_{t}))
#' $$
#' 
#' where $\epsilon$ is the learning rate, $r_t \in \{0, 1, 5\}$ is the size of reinforcement in pennies, and $\rho$ modifies the value of reinforcement. 
#' The intial value for all action weights was set to 0 (or chance responding).
#' 
#' The parameters were estimated using Stan (Stan Development team, 2016) using the `hBayesDM` model with only 
#' [minor modifications](https://github.com/jflournoy/split_bayes/blob/7a9db9add05aca8974aa55bfa8ce7c4231d0f88f/gng_m1_reg.stan).
#' Sampling occured over 1000 iterations in 6 chains, yeilding 3000 samples after discarding the first 50% as the warm-up period. Inspection of the
#' traceplots for the three parameters of interest (below) shows good mixing, though there seems to be a potential issue with chain 3 betwen iterations 700-800.
#'
#' >Guitart-Masip, M., Huys, Q. J. M., Fuentemilla, L., Dayan, P., Duzel, E., & Dolan, R. J. (2012). Go and no-go learning in reward and punishment: Interactions between affect and effect. Neuroimage, 62(1), 154–166.
#' >Stan Development Team. 2016. RStan: the R interface to Stan. R package version 2.14.1.   http://mc-stan.org
#' >Stan Development Team. 2016. Stan Modeling Language Users Guide and Reference Manual, Version 2.15.0.   http://mc-stan.org
#'
#+run models,echo=F,warning=F,message=F,error=T----
if(file.exists(htRW_m1_fname)){
  htRW_m1_fit <- readRDS(htRW_m1_fname)
} else {
  stop('HT chains not found')
}

if(file.exists(dlRW_m1_fname)){
  dlRW_m1_fit <- readRDS(dlRW_m1_fname)
} else {
  stop('DL chains not found')
}

if(file.exists(puRW_m1_fname)){
  puRW_m1_fit <- readRDS(puRW_m1_fname)
} else {
  stop('PU chains not found')
}

stan_trace(htRW_m1_fit, pars = c('mu_p'))
# stan_trace(htRW_m1_fit, pars = c('mu_p','sigma'))
# stan_plot(htRW_m1_fit, pars = c('mu_p'))
# stan_plot(htRW_m1_fit, pars = c('mu_xi', 'mu_ep', 'mu_rho'))
# stan_plot(htRW_m1_fit, pars = c('mu_ep'))
# stan_plot(htRW_m1_fit, pars = c('sigma'))
# stan_plot(htRW_m1_fit, pars = c('ep_pr'))
# stan_plot(htRW_m1_fit, pars = c('rho_pr'))
# stan_plot(htRW_m1_fit, pars = c('xi_pr'))

#'
#' # Mean parameter estimates
#' 
#' ## Epsilon
#' 
#' These are the estimated intervals for the learning rate parameter $\epsilon$.
#' 
#+means1,fig.width=8,fig.height=4
multiplot(stan_plot(htRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Hungry/Thirsty'),
          stan_plot(dlRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Dating/Looking'),
          stan_plot(puRW_m1_fit, pars = c('mu_ep'))+lims(x=c(0,.15))+labs(subtitle='Popular/Unpopular'),
          cols=1)
#' 
#' ## Rho
#' 
#' These are the estimated intervals for the temperature parameter $\rho$.
#' 
#+means2,fig.width=8,fig.height=4
multiplot(stan_plot(htRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Hungry/Thirsty')+lims(x=c(.5,3)),
          stan_plot(dlRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Dating/Looking')+lims(x=c(.5,3)),
          stan_plot(puRW_m1_fit, pars = c('mu_rho'))+labs(subtitle='Popular/Unpopular')+lims(x=c(.5,3)),
          cols=1)
#' 
#' ## Xi
#' 
#' These are the estimated intervals for the noise parameter $\xi$.
#' 
#+means3,fig.width=8,fig.height=4
multiplot(stan_plot(htRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Hungry/Thirsty')+lims(x=c(0,.4)),
          stan_plot(dlRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Dating/Looking')+lims(x=c(0,.4)),
          stan_plot(puRW_m1_fit, pars = c('mu_xi'))+labs(subtitle='Popular/Unpopular')+lims(x=c(0,.4)),
          cols=1)

#'
#'
#+extract_data,echo=F
htRW_m1_fit_extract <- extract(htRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)
dlRW_m1_fit_extract <- extract(dlRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)
puRW_m1_fit_extract <- extract(puRW_m1_fit, pars = c('Qgo', 'Qnogo', 'Wgo', 'Wnogo', 'lp__'), include = FALSE)

rm(htRW_m1_fit, dlRW_m1_fit, puRW_m1_fit);gc()

#'
#'
#+check rho ep cor,echo=F,eval=F
samples <- dim(htRW_m1_fit_extract$ep)[1]
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
#'
#' # Correlations with age
#' 
#' First, make a function to create plots....\
#' 
#+create plot function----
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
#'
#' ## Age cor with ep
#'
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
#'
#' ## Age cor with rho
#'
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

#'
#' ## Age cor with xi
#'
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

#'
#' ## IQ cor with ep
#'
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
#'
#' ## IQ cor with rho
#'
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

#'
#' ## IQ cor with xi
#'
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
#'
#' # Regress parameters on age
#'
#' Helper functions....
#+Age regressions----
library(wesanderson)

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

#'
#' ## Hungry/Thristy parameters
#'
#+Age HT regressions,fig.width=10,fig.height=6----
ep_points_ht <- data_frame(points = apply(htRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2_ht <- regression_predictions(htRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2_ht, points_frame = ep_points_ht, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))

rho_points_ht <- data_frame(points = apply(htRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2_ht <- regression_predictions(htRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2_ht, points_frame = rho_points_ht, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))

xi_points_ht <- data_frame(points = apply(htRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2_ht <- regression_predictions(htRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2_ht, points_frame = xi_points_ht, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))

#'
#' ## Dating/Looking parameters
#'
#+Age DL regressions,fig.width=10,fig.height=6----
ep_points_dl <- data_frame(points = apply(dlRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2_dl <- regression_predictions(dlRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2_dl, points_frame = ep_points_dl, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))

rho_points_dl <- data_frame(points = apply(dlRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2_dl <- regression_predictions(dlRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2_dl, points_frame = rho_points_dl, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))

xi_points_dl <- data_frame(points = apply(dlRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2_dl <- regression_predictions(dlRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2_dl, points_frame = xi_points_dl, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))
#'
#' ## Popular/Unpopular parameters
#'
#+Age PU regressions,fig.width=10,fig.height=6----
ep_points_pu <- data_frame(points = apply(puRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2_pu <- regression_predictions(puRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2_pu, points_frame = ep_points_pu, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))

rho_points_pu <- data_frame(points = apply(puRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2_pu <- regression_predictions(puRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2_pu, points_frame = rho_points_pu, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))

xi_points_pu <- data_frame(points = apply(puRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2_pu <- regression_predictions(puRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2_pu, points_frame = xi_points_pu, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))
#'
#' ## Popular - Hungry contrasted regression
#'
#+Age PU-HT regressions,fig.width=10,fig.height=6----
ep_points_pu_ht <- data_frame(points = apply(puRW_m1_fit_extract$ep_pr - htRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2_pu_ht <- regression_predictions(puRW_m1_fit_extract$ep_pr - htRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2_pu_ht, points_frame = ep_points_pu_ht, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))
rho_points_pu_ht <- data_frame(points = apply(puRW_m1_fit_extract$rho_pr - htRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2_pu_ht <- regression_predictions(puRW_m1_fit_extract$rho_pr - htRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2_pu_ht, points_frame = rho_points_pu_ht, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))
xi_points_pu_ht <- data_frame(points = apply(puRW_m1_fit_extract$xi_pr - htRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2_pu_ht <- regression_predictions(puRW_m1_fit_extract$xi_pr - htRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2_pu_ht, points_frame = xi_points_pu_ht, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))
#'
#' ## Dating - Hungry contrasted regression
#'
#+Age DL-HT regressions,fig.width=10,fig.height=6----
ep_points_dl_ht <- data_frame(points = apply(dlRW_m1_fit_extract$ep_pr - htRW_m1_fit_extract$ep_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
ep_age_poly2_dl_ht <- regression_predictions(dlRW_m1_fit_extract$ep_pr - htRW_m1_fit_extract$ep_pr, stanAges$age)
reg_pred_plot(ep_age_poly2_dl_ht, points_frame = ep_points_dl_ht, palette = regpalette,
              labslist = list(title = "Regressing ep_pr on age across samples", x = 'Age', y = 'ep_pr'))
rho_points_dl_ht <- data_frame(points = apply(dlRW_m1_fit_extract$rho_pr - htRW_m1_fit_extract$rho_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
rho_age_poly2_dl_ht <- regression_predictions(dlRW_m1_fit_extract$rho_pr - htRW_m1_fit_extract$rho_pr, stanAges$age)
reg_pred_plot(rho_age_poly2_dl_ht, points_frame = rho_points_dl_ht, palette = regpalette,
              labslist = list(title = "Regressing rho_pr on age across samples", x = 'Age', y = 'rho_pr'))
xi_points_dl_ht <- data_frame(points = apply(dlRW_m1_fit_extract$xi_pr - htRW_m1_fit_extract$xi_pr, 2, quantile, probs = c(.5)), regressor = stanAges$age)
xi_age_poly2_dl_ht <- regression_predictions(dlRW_m1_fit_extract$xi_pr - htRW_m1_fit_extract$xi_pr, stanAges$age)
reg_pred_plot(xi_age_poly2_dl_ht, points_frame = xi_points_dl_ht, palette = regpalette,
              labslist = list(title = "Regressing xi_pr on age across samples", x = 'Age', y = 'xi_pr'))

#'
#' # Learning through the task
#' 
#+compile learning data----
ht_pr <- as_data_frame(apply(htRW_m1_fit_extract$Pr, c(2,3), median)) %>%
  mutate(id = unique(condition_sum$id)) %>%
  gather(condition_index, Pr,-id) %>%
  mutate(condition_index = as.numeric(sub('V','',condition_index)),
         condition = 'HngT') %>%
  left_join(as_data_frame(apply(htRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .025)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr025,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'HngT')) %>%
  left_join(as_data_frame(apply(htRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .975)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr975,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'HngT'))

dl_pr <- as_data_frame(apply(dlRW_m1_fit_extract$Pr, c(2,3), median)) %>%
  mutate(id = unique(condition_sum$id)) %>%
  gather(condition_index, Pr,-id) %>%
  mutate(condition_index = as.numeric(sub('V','',condition_index)),
         condition = 'DtnL') %>%
  left_join(as_data_frame(apply(dlRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .025)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr025,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'DtnL')) %>%
  left_join(as_data_frame(apply(dlRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .975)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr975,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'DtnL'))

pu_pr <- as_data_frame(apply(puRW_m1_fit_extract$Pr, c(2,3), median)) %>%
  mutate(id = unique(condition_sum$id)) %>%
  gather(condition_index, Pr,-id) %>%
  mutate(condition_index = as.numeric(sub('V','',condition_index)),
         condition = 'PplU') %>%
  left_join(as_data_frame(apply(puRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .025)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr025,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'PplU')) %>%
  left_join(as_data_frame(apply(puRW_m1_fit_extract$Pr, c(2,3), quantile, probs = .975)) %>%
              mutate(id = unique(condition_sum$id)) %>%
              gather(condition_index, Pr975,-id) %>%
              mutate(condition_index = as.numeric(sub('V','',condition_index)),
                     condition = 'PplU'))
prDF <- bind_rows(ht_pr, dl_pr, pu_pr)
prTrialInfoDF <- left_join(prDF, stanSplitDF) %>%
  mutate(condition = factor(condition, levels = c('HngT', 'DtnL', 'PplU'), labels = c('Hngry/Thrsty', 'Dtng/Lkng', 'Pplr/Unpplr')))

id_order_ep <- unique(stanSplitDF$id)[order((ep_points_pu$points+ep_points_dl$points+ep_points_ht$points)/3)]
id_order_rho <- unique(stanSplitDF$id)[order((rho_points_pu$points+rho_points_dl$points+rho_points_ht$points)/3)]

#'
#' ## Individual differences in learning (model estimated)
#'
#' ### Ordered by average learning rate
#'
#+indiv diffs plot,fig.width=18,fig.height=12----
prTrialInfoDF %>%
  filter(!is.na(proportion)) %>%
  mutate(id = factor(id, levels = id_order_ep)) %>%
  group_by(id, condition, proportion) %>%
  mutate(prop_index = 1:n()) %>%
  ggplot(aes(x = trial_index, y = Pr, group = interaction(id, proportion, condition), color = condition)) +
  geom_errorbar(aes(ymin = Pr025, ymax = Pr975), width = 0, size = 1, alpha = .1)+
  geom_line(alpha = .1)+
  geom_point(alpha = .5, size = .5)+
  # geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, bs = 'cr', k = 4), alpha = .15)+
  # geom_line(aes(group = interaction(condition, proportion)), stat = 'smooth', method = 'gam', formula = y ~ s(x, bs = 'cs', k = 4), size = 2)+
  scale_color_manual(name = 'Optimal choice', #limits = c('HngT','DtnL','PplU'),
                     labels = c('Hngry/Thrsty', 'Dtng/Lkng', 'Pplr/Unpplr'),
                     values = wes_palette('Darjeeling')[c(2,1,4)])+
  theme(panel.background = element_blank())+
  facet_wrap(~id)+
  theme(strip.text = element_blank())+
  labs(x = 'Cue-trial number', y = 'Probability of pressing right-arrow')
#'
#' ## Median optimal-choice probability
#'
#' To get a better sense for how the subtle differences in RW parameters are playing out, we can look at this plot of the median probability
#' of choosing the optimal option. This graph collapses across both cues in a domain, and all participants.
#'
#+average prob pressing right----
propDFs <- stanSplitDF %>%
  select(id,proportion,condition_index) %>%
  group_by(condition) %>%
  do({
    data_frame( condDF = list( select(spread(ungroup(.), condition_index, proportion, fill = 1), -id, -condition) ) )
  })

separate_param_cueclass <- function(fit_extract_param, cueclassDF, class_matches){
  require(data.table)
  #in this function we assume that margin 1 indexes the MCMC samples
  indx1 <- which(as.vector(cueclassDF == class_matches[1]), arr.ind = F)
  indx2 <- which(as.vector(cueclassDF == class_matches[2]), arr.ind = F)
  indx1_len <- length(indx1)
  indx2_len <- length(indx2)
  
  nrow_in_mat <- dim(cueclassDF)[1]
  ncol_in_mat <- dim(cueclassDF)[2]
  
  row_id_mat <- matrix(rep(1:nrow_in_mat, ncol_in_mat), nrow=nrow_in_mat)
  
  nsamples <- dim(fit_extract_param)[1]
  # print(nsamples)
  separatedDF <- as.data.table(bind_rows(apply(fit_extract_param, 1, function(x){
    data_frame(Pr = c(as.vector(x)[indx1],as.vector(x)[indx2]),
               cueclass = c(rep(class_matches[1], indx1_len), rep(class_matches[2], indx2_len)))
  })))
   
  separatedDF[, c('row_id', 'sample_id') := list(rep(c(as.vector(row_id_mat)[indx1], as.vector(row_id_mat)[indx2]), nsamples),
                                                 sample_id = rep(1:nsamples, each = (indx1_len+indx2_len)))]
  separatedDF[, cueclass_idx := 1:.N, by = .(cueclass, row_id, sample_id)]
  return(separatedDF)
}

class_matches <- c('80_20','20_80')
ht_pr_trials <- separate_param_cueclass(htRW_m1_fit_extract$Pr, propDFs$condDF[[1]], class_matches)
dl_pr_trials <- separate_param_cueclass(dlRW_m1_fit_extract$Pr, propDFs$condDF[[2]], class_matches)
pu_pr_trials <- separate_param_cueclass(puRW_m1_fit_extract$Pr, propDFs$condDF[[3]], class_matches)

ht_pr_trials[, Pr_across := ifelse(cueclass == '80_20', 1-Pr, Pr)]
dl_pr_trials[, Pr_across := ifelse(cueclass == '80_20', 1-Pr, Pr)]
pu_pr_trials[, Pr_across := ifelse(cueclass == '80_20', 1-Pr, Pr)]

pr_trials <- bind_rows('HngT' = ht_pr_trials, 'DtnL' = dl_pr_trials, 'PplU' = pu_pr_trials, .id = 'condition')
pr_trials_cue_sum <- pr_trials[, as.list(quantile(Pr, probs = c(.025, .5, .975))), by = .(condition, cueclass, cueclass_idx)]
pr_trials_across_cue_sum <- pr_trials[, as.list(quantile(Pr_across, probs = c(.025, .25, .5, .75, .975))), by = .(condition, cueclass_idx)]
rm(pr_trials);gc()
#'
#'
#+median pr plot,fig.width = 6,fig.height=5
ggplot(pr_trials_across_cue_sum, aes(x = cueclass_idx, y = `50%`, group = condition)) +
  geom_line(aes(color = condition), alpha = 1)+
  scale_color_manual(name = 'Condition', limits = c('HngT','DtnL','PplU'),
                     labels = c('Hngry/Thrsty', 'Dtng/Lkng', 'Pplr/Unpplr'),
                     values = wes_palette('Darjeeling')[c(2,1,4)])+
  theme(panel.background = element_blank())+
  labs(x = 'Cue-trial number', y = 'Probability of optimal choice')
