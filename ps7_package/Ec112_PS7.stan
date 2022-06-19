data {
  int<lower=0> N;
  vector[N] log_Mrate;
  vector[N] log_BodySize;
  vector[N] Instar;
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}
model {
  log_Mrate ~ normal(beta0 + beta1 * log_BodySize + beta2 * Instar, sigma);
  beta0 ~ normal(0, 400);
  beta1 ~ normal(0, 400);
  beta2 ~ normal(0, 400);
  sigma ~ gamma(1, 1);
}
generated quantities {
  vector[N] yPred;
  for (i in 1:N) {
    yPred[i] = normal_rng(beta0 + beta1 * log_BodySize[i] + beta2 * Instar[i], sigma);
  }
}

