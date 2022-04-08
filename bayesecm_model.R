library(rstan)

# data includes:
# T, number of time points
# N, number of regions
# n_neighbours, number of neighbours for each region as a vector
# M, maximum number of neighbours for any region
# neighbours, matrix of neighbours: row represents a region, each cell holds the id number of the neighbour
# prices, matrix of rye log-prices: row represents a region, column a date
# missing, matrix indicating missingness: row represents a region, column a date; 0 = present, 1 = missing
par <- readRDS("par.rds")

a <- proc.time()

# initial values
set.seed(1)
k <- sum(par$n_neighbours)
inits <- replicate(4,
                   list(
                     sigma_y = runif(1, 0.03, 0.05),
                     sigma_mu = runif(1, 0.03, 0.05),
                     sigma_a = runif(1, 0.01, 0.03),
                     sigma_c = runif(1, 0.15, 0.25),
                     a = rep(0, par$T),
                     b = abs(rnorm(k, 0.15, 0.1)), 
                     c_std = matrix(0, nrow = k, ncol = par$T),
                     lambda_raw = rep(0, par$N - 1),
                     sigma_lambda = 0.1,
                     phi = 0.6,
                     const_a = 0,
                     mu = par$prices), simplify = FALSE)

model_code <- stan_model("bayesecm_model.stan")

# fit the model
fit <- sampling(model_code, data = par, iter = 8000, warmup = 3000,
                chains = 4, include = FALSE, pars = c("lambda_raw", "c_std"), cores = 4,
                refresh = 10, init = inits, save_warmup = FALSE)

# some results as an example
print(fit, pars = c("sigma_a", "sigma_c", "sigma_mu", "sigma_y"), use_cache = FALSE)

plot(fit, pars = "sigma_y", plotfun = "trace")
plot(fit, pars = "sigma_y", show_density = T)
pairs(fit, pars = c("sigma_mu", "sigma_y"))
