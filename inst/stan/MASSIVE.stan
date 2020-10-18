data {
  int N; // number of observations
  int J; // number of IVs
  int n; // number of Bernoulli trials
  
  // Sufficient statistics (First and second order moments)
  matrix[J + 3, J + 3] SS; 
  
  // Fixed variance for G
  vector<lower = 0>[J] sigma_G;
  
  real mu_X;
  real mu_Y;
  
  // std. dev. for Gaussian prior
  vector<lower = 0>[J] sigma_gamma;
  vector<lower = 0>[J] sigma_alpha;
}


parameters {
  
  real sbeta;
  vector[J] sgamma;
  vector[J] salpha;

  real<lower = 0> skappa_X;
  real skappa_Y;
  
  real log_sigma_X;
  real log_sigma_Y;
}

transformed parameters {
  
  // real skappa_X = fabs(ur_skappa_X);
  real sigma_X = exp(log_sigma_X);
  real sigma_Y = exp(log_sigma_Y);
  
  // original parameters (before rescaling)
  real beta = sbeta * sigma_Y / sigma_X;
  vector[J] gamma = sgamma * sigma_X ./ sigma_G;
  vector[J] alpha = salpha * sigma_Y ./ sigma_G;
  
  real kappa_X = skappa_X * sigma_X;
  real kappa_Y = skappa_Y * sigma_Y;
  real prod_skappa = skappa_X * skappa_Y;
  
  
  /* Covariance of conditional Gaussian distribution */
  matrix[2, 2] sigma;
  matrix[J + 3, 2] f;
  
  sigma[1, 1] = sigma_X^2 * (1 + skappa_X^2);
  sigma[1, 2] = sigma_X * sigma_Y * (sbeta * (1 + skappa_X^2) + skappa_X * skappa_Y);
  sigma[2, 1] = sigma[1, 2];
  sigma[2, 2] = sigma_Y^2 * (1 + sbeta^2 + (skappa_Y + sbeta * skappa_X)^2);


  f[1, 1] = - mu_X;
  f[1, 2] = - (mu_Y + beta * mu_X);
  
  for (j in 1:J) {
    f[j + 1, 1] = - sigma_X * sgamma[j] / sigma_G[j];
    f[j + 1, 2] = - sigma_Y * (salpha[j] + sbeta * sgamma[j]) / sigma_G[j];
  }
  
  f[J + 2, 1] = 1;
  f[J + 2, 2] = 0;
  f[J + 3, 1] = 0;
  f[J + 3, 2] = 1;
}


model {
  
  // NOTE: in order to use bridge sampling, we need to preserve constants
  
  /* Prior on Structural Parameters */
  target += normal_lpdf(sgamma | 0, sigma_gamma);
  target += normal_lpdf(salpha | 0, sigma_alpha);
  target += normal_lpdf(sbeta | 0, 10);

  /* Prior on Confounding Coefficients */
  target += normal_lpdf(skappa_X | 0, 10) + log(2); // half_normal
  target += normal_lpdf(skappa_Y | 0, 10);
  
  // conditional likelihood, not full likelihood, i.e. p(X, Y | G)
  target += N * (- 0.5 * (log(4 * pi() * pi()) + log_determinant(sigma) + 
  trace(inverse(sigma) * f' * SS * f)));
  
}

