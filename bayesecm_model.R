library(rstan)

# data includes:
# T, number of time points
# N, number of regions
# n_neighbours, number of neighbours for each region as a vector
# M, maximum number of neighbours for one region
# neighbours, matrix of neighbours
# prices, matrix of rye log-prices
# missing, matrix indicating missingness
par <- readRDS("par.rds")

# initial values
set.seed(1)
k <- sum(par$n_neighbours)
inits <- replicate(4,
                   list(
                     sigma_y = runif(1, 0.05, 0.5),
                     sigma_mu = runif(1, 0.1, 1),
                     sigma_a = runif(1, 0.01, 0.1),
                     sigma_c = runif(1, 0.1, 1),
                     a = rep(0, par$T),
                     b = abs(rnorm(k, 0.15, 0.1)),
                     c_std = matrix(rnorm(k * par$T, -2, 0.5), nrow = k, ncol = par$T),
                     mu = par$prices), simplify = FALSE)

model_code <- stan_model("bayesecm_model.stan")

# fit the model
fit <- sampling(model_code, data = par, iter = 6000, warmup = 3000,
                     chains = 4, include = FALSE, pars = "c_std", cores = 4,
                     refresh = 50, init = inits, save_warmup = FALSE)

# some results as an example
print(fit, pars = c("sigma_a", "sigma_c", "sigma_mu", "sigma_y"))

plot(fit, pars = "sigma_y", plotfun = "trace")
plot(fit, pars = "sigma_y", show_density = T)
pairs(fit, pars = c("sigma_mu", "sigma_y"))
