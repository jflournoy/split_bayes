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
}'

N=5

somedata <- list(N=N,
                 x=runif(N,0,1))

fit <- stan(model_code = model, data = somedata, algorithm = 'Fixed_param',  iter = 1000, chains = 8)

