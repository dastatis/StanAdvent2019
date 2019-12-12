rm(list=ls())
# Functions ---------------------------------------------------------------

## inv_logit
ilogit <- function(x) {
  re <- 1/(1+exp(-x))
  return(re)
}

## log_sum_exp
lse <- function(x) {
  re <- exp(sum(log(x)))
  return(re)
}

## ggplot ( scatter plot )
ggpoint <- function(x,y) {
  re <- 
    data.frame(x = x, y = y) %>% 
    ggplot(aes(x = .[,1], y = .[,2])) +
    geom_point() +
    geom_abline(slope=1) +
    xlab("") +
    ylab("") +
    theme_bw()
  return(re)
}

# Packages ----------------------------------------------------------------

library(tidyverse)
library(extraDistr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# source("R/GeneralRFunc.R")

# Data and Prior settings -------------------------------------------------

I <- 1000 ##  # of Respondents
J <- 10 ## # of Items

set.seed(3456)
alpha <- rlnorm(J, meanlog = 0, sdlog = 1) ## discrimination
# alpha <- extraDistr::rhcauchy(J, sigma = 2.5) ## discrimination
# alpha <- stats::runif(J, min = 0.5, max = 2.5) ## discrimination
beta <- stats::rnorm(J, mean = 0, sd = 2) ## difficulty
gamma <- stats::runif(J, min = 0, max = 0.2) ## guessing 
upsilon <- stats::runif(J, min = 0.8, max = 1) ## 
zeta <- extraDistr::rhnorm(J, sigma = 1) ## asymmetric
# zeta <- qnorm(runif(J,pnorm(0, mean = 1, sd = 1, lower.tail = 1), pnorm(Inf, mean = 1, sd = 1)), mean = 1, sd = 1) ## asymmetric
theta <- rnorm(I, mean = 0, sd = 1) ## trait or ability

# Sampling settings --------------------------------------------------------

chain <- 4
iter <- 500
warmup <- floor(iter/2)
thin <- 1
seed <- 1234
output_samples <- chain * ((iter-warmup)/thin)

## watanabe beta (invese temparture) setting
mode = 0  ## 0 = usual MCMC, 1 = WBIC MCMC

# 1PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- ilogit( theta[i]-beta[j] )
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])
  }
}

## thetaの各値で反応確率を計算(作図やIIF,TIFを算出するため)
pseudo_theta <- seq(from = -10, to = 10, length.out =  10000)
pseudo_prob <- matrix(NA, nrow = 10000, ncol = J)
for (i in 1:10000) {
  for (j in 1:J) {
    pseudo_prob[i,j] <- ilogit( pseudo_theta[i]-beta[j] )
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_1pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("~/Documents/git/BayesianPsychologicalModeling/Model/IRT/IRT_1PL.stan")

fit_1pl <- sampling(
  object = model,
  data = standat_1pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_1pl)

est_beta <- get_posterior_mean(fit_1pl,"beta")[,"mean-all chains"]

rc_p_beta <- ggpoint(beta,est_beta)
rc_p_beta

# 2PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- ilogit( alpha[j] * (theta[i]-beta[j]) )
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])  }
}

## thetaの各値で反応確率を計算(作図やIIF,TIFを算出するため)
pseudo_N <- 1000
pseudo_theta <- seq(from = -10, to = 10, length.out =  pseudo_N)
pseudo_prob <- matrix(NA, nrow = pseudo_N, ncol = J)
for (i in 1:pseudo_N) {
  for (j in 1:J) {
    pseudo_prob[i,j] <- ilogit( alpha[j] * (pseudo_theta[i]-beta[j]) )
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_2pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("Model/IRT/IRT_2PL.stan")

fit_2pl <- sampling(
  object = model,
  data = standat_2pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_2pl)

est_alpha <- get_posterior_mean(fit_2pl,"alpha")[,"mean-all chains"]
est_beta <- get_posterior_mean(fit_2pl,"beta")[,"mean-all chains"]

rc_p_alpha <- ggpoint(alpha,est_alpha)
rc_p_beta <- ggpoint(beta,est_beta)

# 3PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- gamma[j] + ( (1 - gamma[j]) * ilogit( alpha[j] * (theta[i]-beta[j]) ) )
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])
  }
}

## thetaの各値で反応確率を計算(作図やIIF,TIFを算出するため)
pseudo_N <- 1000
pseudo_theta <- seq(from = -10, to = 10, length.out =  pseudo_N)
pseudo_prob <- matrix(NA, nrow = pseudo_N, ncol = J)
for (i in 1:pseudo_N) {
  for (j in 1:J) {
    pseudo_prob[i,j] <- gamma[j] + ( (1 - gamma[j]) * ilogit( alpha[j] * (pseudo_theta[i]-beta[j]) ) )
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_3pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("Model/IRT/IRT_3PL.stan")

fit_3pl <- sampling(
  object = model,
  data = standat_3pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_3pl)

est_alpha <- get_posterior_mean(fit_3pl,"alpha")[,"mean-all chains"]
est_beta <- get_posterior_mean(fit_3pl,"beta")[,"mean-all chains"]
est_gamma <- get_posterior_mean(fit_3pl,"gamma")[,"mean-all chains"]

rc_p_alpha <- ggpoint(alpha,est_alpha)
rc_p_beta <- ggpoint(beta,est_beta)
rc_p_gamma <- ggpoint(gamma,est_gamma)

# 4PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- gamma[j] + ( (upsilon[j] - gamma[j]) * ilogit( alpha[j] * (theta[i]-beta[j]) ) )
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])
  }
}

## thetaの各値で反応確率を計算(作図やIIF,TIFを算出するため)
pseudo_N <- 1000
pseudo_theta <- seq(from = -10, to = 10, length.out =  pseudo_N)
pseudo_prob <- matrix(NA, nrow = pseudo_N, ncol = J)
for (i in 1:pseudo_N) {
  for (j in 1:J) {
    pseudo_prob[i,j] <- gamma[j] + ( (upsilon[j] - gamma[j]) * ilogit( alpha[j] * (pseudo_theta[i]-beta[j]) ) )
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_4pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("~/Documents/git/BayesianPsychologicalModeling/Model/IRT/IRT_4PL.stan")

fit_4pl <- sampling(
  object = model,
  data = standat_4pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_4pl)

est_alpha <- get_posterior_mean(fit_4pl,"alpha")[,"mean-all chains"]
est_beta <- get_posterior_mean(fit_4pl,"beta")[,"mean-all chains"]
est_gamma <- get_posterior_mean(fit_4pl,"gamma")[,"mean-all chains"]
est_upsilon <- get_posterior_mean(fit_4pl,"upsilon")[,"mean-all chains"]

rc_p_alpha <- ggpoint(alpha,est_alpha)
rc_p_beta <- ggpoint(beta,est_beta)
rc_p_gamma <- ggpoint(gamma,est_gamma)
rc_p_upsilon <- ggpoint(upsilon,est_upsilon)

rc_p_alpha
rc_p_beta
rc_p_gamma
rc_p_upsilon

# 5PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- gamma[j] + ( (upsilon[j] - gamma[j]) * (ilogit( alpha[j] * (theta[i]-beta[j]) )^zeta[j] ) )
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])
    
  }
}

## thetaの各値で反応確率を計算(作図やIIF,TIFを算出するため)
pseudo_N <- 1000
pseudo_theta <- seq(from = -10, to = 10, length.out =  pseudo_N)
pseudo_prob <- matrix(NA, nrow = pseudo_N, ncol = J)
for (i in 1:pseudo_N) {
  for (j in 1:J) {
    pseudo_prob[i,j] <- gamma[j] + ( (upsilon[j] - gamma[j]) * (ilogit( alpha[j] * (pseudo_theta[i]-beta[j]) ))^zeta[j] )
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_5pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("IRT_5PL.stan")

fit_5pl <- sampling(
  object = model,
  data = standat_5pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_5pl)

est_alpha <- get_posterior_mean(fit_5pl,"alpha")[,"mean-all chains"]
est_beta <- get_posterior_mean(fit_5pl,"beta")[,"mean-all chains"]
est_gamma <- get_posterior_mean(fit_5pl,"gamma")[,"mean-all chains"]
est_upsilon <- get_posterior_mean(fit_5pl,"upsilon")[,"mean-all chains"]
est_zeta <- get_posterior_mean(fit_5pl,"zeta")[,"mean-all chains"]

rc_p_alpha <- ggpoint(alpha,est_alpha)
rc_p_beta <- ggpoint(beta,est_beta)
rc_p_gamma <- ggpoint(gamma,est_gamma)
rc_p_upsilon <- ggpoint(upsilon,est_upsilon)
rc_p_zeta <- ggpoint(zeta,est_zeta)

rc_p_alpha
rc_p_beta
rc_p_gamma
rc_p_upsilon
rc_p_zeta



# Asymmetric 2PL IRT -----------------------------------------------------------------
prob <- matrix(NA, nrow = I, ncol = J)
y <- matrix(NA, nrow = I, ncol = J)
pseudo_theta <- seq(from = -10, to = 10, length.out =  I)
pseudo_prob <- matrix(NA, nrow = I, ncol = J)
for (i in 1:I) {
  for (j in 1:J) {
    prob[i,j] <- ilogit( alpha[j] * (theta[i]-beta[j]) )^zeta[j]
    pseudo_prob[i,j] <- (ilogit( alpha[j] * (pseudo_theta[i]-beta[j]) ))^zeta[j]
    y[i,j] <- extraDistr::rbern(n = 1, prob = prob[i,j])
  }
}

## Data Visualizations
# Item Characteristic Curve
pseudo_prob %>% 
  data.frame(prob = .) %>% 
  bind_cols(theta = pseudo_theta) %>% 
  gather(key = "item", value = "prob", -"theta") %>% 
  mutate(
    item = as.integer(str_replace(item, "^prob.", ""))
  ) %>% 
  filter(item <= 10) %>% 
  ggplot(aes(x = theta, y = prob, color = factor(item))) +
  geom_line(size = 1) +
  theme_minimal()

standat_Asym2pl <- list(
  I = I,
  J = J,
  y = y,
  mode = mode
)

model <- stan_model("Model/IRT/IRT_Asymmetric2PL.stan")

fit_Asym2pl <- sampling(
  object = model,
  data = standat_Asym2pl,
  chains = chain,
  iter = iter,
  warmup = warmup,
  thin = thin,
  seed = seed
)

stan_convergence(fit_Asym2pl)

est_alpha <- get_posterior_mean(fit_Asym2pl,"alpha")[,"mean-all chains"]
est_beta <- get_posterior_mean(fit_Asym2pl,"beta")[,"mean-all chains"]
est_zeta <- get_posterior_mean(fit_Asym2pl,"zeta")[,"mean-all chains"]

est_alpha <- get_posterior_map(fit_Asym2pl,"alpha")
est_beta <- get_posterior_map(fit_Asym2pl,"beta")
est_zeta <- get_posterior_map(fit_Asym2pl,"zeta")

rc_p_alpha <- ggpoint(alpha,est_alpha)
rc_p_beta <- ggpoint(beta,est_beta)
rc_p_zeta <- ggpoint(zeta,est_zeta)

rc_p_alpha
rc_p_beta
rc_p_zeta

