# A Bayesian spatio-temporal error correction analysis of markets during the Finnish 1860s famine

Materials for the paper _A Bayesian spatio-temporal error correction analysis of markets during the Finnish 1860s famine_ by Tiia-Maria Pasanen, Miikka Voutilainen, Jouni Helske and Harri Högmander

The model is defined in the file ```bayesecm_model.stan``` and used in the corresponding R script found from ```bayesecm_model.R```. The data used is in the file ```par.rds```, which is read into R in the file ```bayesecm_model.R```. The data includes the monthly rye prices in log-scale from January 1861 to December 1869, and each row corresponds one region and the rownumber the id number of that region. Also the neighbourhood matrix is included.
