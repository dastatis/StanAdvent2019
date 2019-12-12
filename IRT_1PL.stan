data{
  int<lower = 0> I;
  int<lower = 0> J;
  int<lower = 0, upper = 1> y[I,J];
  int<lower = 0, upper = 1> mode;
}

transformed data {
  real<lower = 0, upper = 1> watanabe_beta;
  if ( mode == 0 ) {
    watanabe_beta = 1.0;
  }
    else { // mode == 1
    watanabe_beta = 1.0/log(I);
  }
}

parameters {
  vector[J] beta;
  vector[I] theta;
}

transformed parameters {
  real<lower = 0, upper = 1> prob[I,J];
  for (i in 1:I) {
    for (j in 1:J) {
      prob[i,j] = inv_logit((theta[i] - beta[j]));
    }
  }
}

model{
  target += normal_lpdf(beta | 0 , 5);
  target += normal_lpdf(theta | 0 , 1);
  for (i in 1:I) {
    for (j in 1:J) {
      target += watanabe_beta * bernoulli_lpmf(y[i,j] | prob[i,j]);
    }
  }
}

generated quantities {
  int<lower = 0, upper = 1> y_pred[I,J];
  vector[I] log_lik;
  vector[J] log_lik_temp;
  for (i in 1:I) {
    for (j in 1:J) {
      log_lik_temp[j] = bernoulli_lpmf(y[i,j] | prob[i,j]);
      y_pred[i,j] = bernoulli_rng(prob[i,j]);
    }
    log_lik[i] = log_sum_exp(log_lik_temp);
  }
}
