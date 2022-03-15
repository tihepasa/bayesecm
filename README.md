# A Bayesian spatio-temporal analysis of markets during the Finnish 1860s famine

Materials for the paper *A Bayesian spatio-temporal analysis of markets during the Finnish 1860s famine* by Tiia-Maria Pasanen, Miikka Voutilainen, Jouni Helske and Harri Högmander

The model is defined in the file `bayesecm_model.stan` and used in the corresponding R script found from `bayesecm_model.R`. The data used is in the file `par.rds`, which is read into R in the file `bayesecm_model.R`. The data includes the monthly regional rye prices on log scale from January 1861 to December 1869, and each row corresponds one region and the rownumber the id number of that region. Also the neighbourhood matrix is included. The overall 80 Finnish regions in our context correspond to 48 administrational districts and 32 town administrations of that time. The regions can be seen from the figures in the folder `figs`, which includes illustrations of regional price increases simulated according to the model.

## Comparative analyses

Our model was compared with a model with different prior distributions and a more classic SAR model, which are described below.

### Effect of priors

The priors we use for our model were chosen to reduce computational burden. The model we present is also fitted with wider priors, which can be seen from the corresponding Stan file `bayesecm_wider_priors.stan`. The priors do not have strong influence on the results. The results are presented in the table and figures below.

### Comparison with classic SAR

A simpler SAR model was also fitted to the same data to compare the latent level outputs $\mu_{i,t}$. In this case we defined $\mu_t - \mu_{t-1} = \rho W (\mu_t - \mu_{t-1}) + \epsilon$, leading to $`\mu_t \sim N(\mu_{t-1},(I -\rho W)^{-1}((I -\rho W)^{-1})^T`$, where $W$ is the fixed adjacency matrix after spectral normalization, and $-1 < \rho < 1$ is the unknown spatial dependency parameter. Thus in addition to latent $\mu$, this model contained only three unknown parameters, $\sigma_y$, $\sigma_{\mu}$, and $\rho$, which were estimated as 0.034 ([0.032, 0.035]), 0.055 ([0.054, 0.057]) and 0.778 ([0.755, 0.801]) respectively, implying strong spatial dependence and larger unexplained variation in the latent log-prices than our main model. The estimates of the $\mu_{i,t}$ (plotted below) resembled those of our main model except that the posterior intervals were wider with reference to the main model.
