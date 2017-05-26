## ---- Load data ----

library(tidyverse)
library(stringr)

excl_ids <- c(110) # 110 didn't play the same game
root_dir <- '/home/flournoy/data/split_bayes' # This would be '/Volumes/' if you've mounted locally
raw_data_dir <- file.path(root_dir, 'experiment_data/')
processed_data_dir <- file.path(root_dir, 'processed/')
#previous iterations of this task included 'BannasOranges' and we don't want those trials
conditions <- c('HungryThirsty', 'DatingLooking', 'PopularUnpopular')
#We will use this later to separate the metacog questions from trial behavior
metacog_condition <- 'Confidence'
file_match_pattern <- 'split.*csv'

age_data <- read_csv(file.path(root_dir, '/age_gender_iq.csv'))

chain_dir <- '/data/jflournoy/split/bayes/'
htRW_m1_fname <- file.path(chain_dir, 'htRW_m1_stan.RDS')
dlRW_m1_fname <- file.path(chain_dir, 'dlRW_m1_stan.RDS')
puRW_m1_fname <- file.path(chain_dir, 'puRW_m1_stan.RDS')

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


message(paste0('Files: ', length(unique(files$file))))

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
  mutate(proportion = condition, 
         condition=paste(sort(unlist(str_split(context, '_'))), collapse=''), 
         condition=ifelse(trial_id=='metacog', 'Confidence', condition),
         correcttrue=ifelse(trial_id=='metacog', 
                            as.numeric(score_response), correcttrue)) %>%
  filter(condition %in% c(conditions, metacog_condition)) %>%
  group_by(id, condition) %>%
  mutate(condition_trial_index=1:n(),
         block=((condition_trial_index-1) %/% 16) +1) %>%
  ungroup() %>%
  select(id, correcttrue, pressed_r, outcome, condition, condition_trial_index, sex, stim_image, trial_index, block, proportion)

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
