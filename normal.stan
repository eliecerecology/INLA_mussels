data {
  int<lower=1> N; //the number of observations
  int<lower=1> K; //number of predictiors
  matrix[N,K] x; //predictor matrix
  vector[N] y; //the response variable
}
parameters {
  real alpha; //intercept         
  vector[K] beta; //matrix of group-level regression coefficients
  real<lower=0> sigma; //standard deviation of the individual observations
}
model {
  y ~ normal(alpha + x*beta, sigma); // likelihood
}
