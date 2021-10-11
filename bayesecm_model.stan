// Stan code to implement the model for the rye prices

functions {

  // custom function for the likelihood
  real likelihood(matrix mu, matrix eta, real sigma, matrix I_B) {
    
    int N = rows(mu);
    int T = cols(mu);
    real inv_sigma2 =  1.0 / square(sigma);
    
    real ll = (T - 1) * (- N * log(sigma) + log_determinant(I_B));
    
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
  
  n_neighbours_cumsum[1] = 1;
  for (i in 2:(N + 1)) n_neighbours_cumsum[i] = n_neighbours_cumsum[i - 1] + n_neighbours[i - 1];
  
}
parameters {
  vector<lower=0>[N] lambda; // site specific trend log-intensity
  vector[T] a; // common trend
  vector<lower=0>[n] b; // short-term effect
  matrix[n, T] c_std; // raw value for the reparameterized long-term effect i.e. error correction coefficient
  
  real const_a; // constant for AR(1) process
  real<lower=0, upper=1> phi; // coefficient for AR(1) process
  
  // standard deviations
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_a;
  real<lower=0> sigma_c; // for the actual (reparameterized) gamma
  real<lower=0> sigma_y; 
  
  matrix[N, T] mu; // latent, unobserved (log-)prices
  real<lower=0> sigma_mu;
  
}
transformed parameters {
  
  matrix[n, T] c; // actual reparameterized error correction coefficients
  real log_p_mu; // log-density of mu
  
  // reparameterize random walk gamma
  c[, 1] = -2 + 2 * c_std[, 1]; // N(-2, 4)
  for (t in 2:T) c[, t] = c[, t - 1] + sigma_c * c_std[, t];
  c = -inv_logit(c); // restrict -1 < gamma < 0
  
  // eta
  {
    matrix[N, T] eta;
    
    // (I - B)
    matrix[N, N] I_B = rep_matrix(0, N, N);
    
    for (i in 1:N) {
      int idx_alku = n_neighbours_cumsum[i];
      int idx_loppu = n_neighbours_cumsum[i+1] - 1;
  
      I_B[neighbours[i, 1:n_neighbours[i]], i] = -(b[idx_alku:idx_loppu]);
      I_B[i, i] += 1;
    }
    
    // D
    for (t in 1:T) {
      matrix[N, N] D = rep_matrix(0, N, N);
      
      for (i in 1:N) {
        int idx_alku = n_neighbours_cumsum[i];
        int idx_loppu = n_neighbours_cumsum[i + 1] - 1;
        
        D[i, i] = 1 + sum(c[idx_alku:idx_loppu, t]); // diagonal
        D[neighbours[i, 1:n_neighbours[i]], i] = -(b[idx_alku:idx_loppu] + c[idx_alku:idx_loppu, t]); // non-zero elements
      }
      
      eta[, t] = a[t] * lambda + D' * mu[, t];
    }

  log_p_mu = likelihood(mu, eta, sigma_mu, I_B');
  }
  
}
  
model {
  
  // actual and latent prices
  sigma_y ~ gamma(2, 20);
  sigma_mu ~ gamma(2, 20);
  mu[, 1] ~ normal(3, 0.5);
  
  // AR(1) alpha
  a[1] ~ normal(const_a / (1 - phi), sigma_a / sqrt(1 - phi^2));
  a[2:T] ~ normal(const_a + phi * a[1:(T - 1)], sigma_a);
  sigma_a ~ gamma(2, 100);
  const_a ~ normal(0, 0.1);
  phi ~ beta(2, 2);
  
  // lambda ~ N(1, sigma_lambda^2)
  sigma_lambda ~ normal(0, 0.5);
  lambda ~ normal(1, sigma_lambda);
  
  // beta
  b ~ gamma(0.5, 2); // 95% prior interval [0.00025, 1.25597], median 0.11
  
  // gamma
  to_vector(c_std) ~ std_normal();
  sigma_c ~ gamma(2, 10);
  
  target += log_p_mu;
  
  for (t in 1:T) {
    for (i in 1:N) {
      if (missing[i, t] == 0) {
        prices[i, t] ~ normal(mu[i, t], sigma_y);
      }
    }
  }
  
}
