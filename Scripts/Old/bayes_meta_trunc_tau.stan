// bayes_meta_trunc_tau.stan
data {
  int<lower=1> K;              // number of studies
  vector[K] y;                 // study-level estimates
  vector<lower=0>[K] sigma;    // known standard errors
}

parameters {
  real mu;                                 // overall mean effect
  real<lower=0.2, upper=0.5> tau;          // heterogeneity SD, truncated
  vector[K] theta;                         // study-specific true effects
}

model {
  // Priors
  mu ~ normal(0, 4);                       // adjust if you like
  // tau has an implicit Uniform(0.2, 0.5) prior via its bounds

  // Hierarchical model
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}

generated quantities {
  real tau2;
  real theta_new;

  tau2 = square(tau);

  // Posterior predictive for a new study's true effect
  theta_new = normal_rng(mu, tau);
}
