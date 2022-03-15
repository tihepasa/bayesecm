// Stan code to implement the model for the rye prices

functions {

  // custom function for the likelihood
  real likelihood(matrix mu, matrix eta, real sigma, matrix I_B, real rho, vector lambda) {
    
    int N = rows(mu);
    int T = cols(mu);
    real inv_sigma2 =  1.0 / square(sigma);
    
    real logdet = sum(log1m(rho * lambda));
    
    real ll = (T - 1) * (- N * log(sigma) + logdet);
    
    for (t in 1:(T - 1)) {
      ll += -0.5 * dot_self(I_B * mu[, t + 1] - eta[, t]) * inv_sigma2;
    }
    return ll;
  }
  
}
data {
  
  int<lower=0> T; // number of time points
  int<lower=0> N; // number of regions
  int<lower=0> n_neighbours[N]; // number of neighbours for each region (cannot be a neighbour of itself)
  int<lower=0> M; // maximum number of neighbours
  int<lower=0> neighbours[N, M]; // matrix of neighbours
  matrix[N, T] prices; // matrix of (log-)prices
  
  matrix[N, T] missing; // matrix indicating missingness

}
transformed data {
  
  // helper variables for looping through the coefficient matrices
  
  int<lower=0> n_neighbours_cumsum[N + 1]; // vector of cumulative sum of regional number of neighbours
  int<lower=0> n = sum(n_neighbours); // total number of neighbour pairs
  matrix[N, N] W = rep_matrix(0, N, N);
  vector[N] lambda;
  
  n_neighbours_cumsum[1] = 1;
  for (i in 2:(N + 1)) n_neighbours_cumsum[i] = n_neighbours_cumsum[i - 1] + n_neighbours[i - 1];
  
//  for (i in 1:N) {
//    W[neighbours[i, 1:n_neighbours[i]], i] += rep_vector(1, n_neighbours[i]);
//    W[, i] /= sum(W[, i]); // each row (column) sums to 1 to ensure invertibility
//  }
//  W = W';
  
  for (i in 1:N) {
    W[neighbours[i, 1:n_neighbours[i]], i] += 1.0;
  }
  W = W / max(fabs(eigenvalues_sym(W)));
  lambda = eigenvalues_sym(W);
}
parameters {
  
//  vector<lower = 0>[N] lambda; // site specific trend log-intensity
//  vector[T] a; // common trend
//  vector<lower=0>[n] b; // short-term effect
  real<lower = -1, upper = 1> rho; // short-term effect
//  matrix[n, T] c_std; // raw value for the reparameterized long-term effect i.e. error correction coefficient
//  vector<lower = -1, upper = 0>[n] c; // raw value for the reparameterized long-term effect i.e. error correction coefficient
  
//  real const_a; // constant for AR(1) process
//  real<lower=0, upper=1> phi; // coefficient for AR(1) process
  
  // standard deviations
//  real<lower=0> sigma_lambda;
//  real<lower=0> sigma_a;
//  real<lower=0> sigma_c; // for the actual (reparameterized) gamma
  real<lower=0> sigma_y; 
  
  matrix[N, T] mu; // latent, unobserved (log-)prices
  real<lower=0> sigma_mu;
  
}
transformed parameters {
  real log_p_mu = likelihood(mu, mu, sigma_mu, add_diag(-rho * W, 1), rho, lambda);
}
model {
  
  // actual and latent prices
  sigma_y ~ gamma(2, 20); // wider to avoid multimodality
  sigma_mu ~ gamma(2, 20);
  mu[, 1] ~ normal(3, 0.5);
  
  // AR(1) alpha
//  a[1] ~ normal(const_a / (1 - phi), sigma_a / sqrt(1 - phi^2));
//  a[2:T] ~ normal(const_a + phi * a[1:(T - 1)], sigma_a);
//  sigma_a ~ gamma(2, 100);
//  const_a ~ normal(0, 0.1);
//  phi ~ beta(2, 2);
  
  // lambda ~ N(1, sigma_lambda^2)
//  sigma_lambda ~ normal(0, 0.5);
//  lambda ~ normal(1, sigma_lambda);
  
  // beta
//  b ~ gamma(0.5, 2); // 95% prior interval [0.00025, 1.25597], median 0.11
// uniform(0, inf) prior for rho i.e. scalar beta
  
  // gamma
//  c ~ beta(2, 2);
//  sigma_c ~ gamma(2, 10);
  
  target += log_p_mu;
  
  for (t in 1:T) {
    for (i in 1:N) {
      if (missing[i, t] == 0) {
        prices[i, t] ~ normal(mu[i, t], sigma_y);
      }
    }
  }
  
}
generated quantities {
  
  matrix[N, T] ypred; // for visual posterior predictive checks etc
  matrix[N, T] log_py_pred; // for model comparisons
  
  for (t in 1:T) {
    for (i in 1:N) {
      // sample from predictive distribution
      ypred[i, t] = normal_rng(mu[i, t], sigma_y);
      // log-density for test data (compute for every element for simplicity)
      log_py_pred[i, t] = normal_lpdf(prices[i, t] | mu[i, t], sigma_y);
    }
  }
  
}
