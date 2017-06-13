#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Associative Learning Simulation"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("xi",
                     "Xi",
                     min = -2,
                     max = 2,
                     value = 0,
                     step=.05),
         sliderInput("ep",
                     "Epsilon",
                     min = -2,
                     max = 2,
                     value = 0,
                     step=.05),
         sliderInput("b",
                     "Bias",
                     min = -2,
                     max = 2,
                     value = 0,
                     step=.05),
         sliderInput("rho",
                     "Rho",
                     min = -2,
                     max = 2,
                     value = 0,
                     step=.05)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("trialsPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  p_right <- data.frame(expand.grid(cue=1:6, reward = c(1,5)), pcorrect_if_pressed_r=c(rep(.2,3), rep(.8,3)))
  
  Trials <- p_right[sample(1:dim(p_right)[1], size = 8*48, replace = T),]
  Trials$crct_if_right <- rbinom(dim(Trials)[1], size = 1, prob = Trials$pcorrect_if_pressed_r)
  Trials$outcome_r <- Trials$crct_if_right*Trials$reward
  Trials$outcome_l <- (1-Trials$crct_if_right)*Trials$reward
  
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
      scale_shape_manual(values = c(25,24), name = 'Last press was...', breaks = c(1,0), labels = c('Right', 'left')) +
      scale_y_continuous(breaks = c(0,1), labels = c('left', 'right'))+
      geom_point() + 
      geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, k = 8,  bs = "cs")) + 
      facet_wrap(pcorrect_if_pressed_r~cue, nrow = 2)+
      theme(panel.background = element_blank(), strip.text = element_blank())
    return(aplot)
  }
  
  output$trialsPlot <- renderPlot({
      plot_RW_run(trials = Trials,
                  mu_p = c(xi = input$xi, ep = input$ep, b = input$b, rho = input$rho))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

