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
  vector<lower = 0>[J] alpha;
  vector[J] beta;
  vector<lower = 0.0, upper = 0.2>[J] gamma;
  vector[I] theta;
}

transformed parameters {
  real<lower = 0, upper = 1> prob[I,J];
  for (i in 1:I) {
    for (j in 1:J) {
      prob[i,j] = gamma[j] + (1-gamma[j]) * inv_logit(alpha[j] * (theta[i] - beta[j]));
    }
  }
}

model{
  target += lognormal_lpdf(alpha | 0 , 1) - lognormal_lccdf(0 | 0, 1);
  target += normal_lpdf(beta | 0 , 5);
  target += uniform_lpdf(gamma | 0.0 , 0.2) - uniform_lccdf(0 | 0.0 , 0.2) -  uniform_lcdf(0.2 | 0.0 , 0.2);
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
