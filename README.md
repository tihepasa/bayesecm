# A Bayesian spatio-temporal analysis of markets during the Finnish 1860s famine

Materials for the paper *A Bayesian spatio-temporal analysis of markets during the Finnish 1860s famine* by Tiia-Maria Pasanen, Miikka Voutilainen, Jouni Helske and Harri Högmander

The paper itself is available at https://arxiv.org/pdf/2106.06268.pdf.

The model is defined in the file `bayesecm_model.stan` and used in the corresponding R script found from `bayesecm_model.R`. The data used is in the file `par.rds`, which is read into R in the file `bayesecm_model.R`. The data includes the monthly regional rye prices on log scale from January 1861 to December 1869, and each row corresponds one region and the rownumber the id number of that region. Also the neighbourhood matrix is included. The overall 80 Finnish regions in our context correspond to 48 administrational districts and 32 town administrations of that time. The regions can be seen from the figures in the folder `figs`, which includes illustrations of regional price increases simulated according to the model.

## Comparative analyses

The material below is intended to represent some results of the comparative analysis we have done to our model. A model with different prior distributions and a more classic SAR model are considered here. In order to see the results of our main model, see our paper.

### Effect of priors

The priors we use for our model were chosen to reduce computational burden. The model we present is also fitted with wider priors, which can be seen from the corresponding Stan file `bayesecm_wider_priors.stan`. The priors do not have strong influence on the results. The results are presented in the table and figures below.

**Table 1.** Scalar parameters of our model with stricter priors and corresponding parameters of the same model with wider priors.
| Parameter                  | Strict priors |             |        |          |           | Wide priors |             |        |          |           |
|----------------------------|---------------|-------------|--------|----------|-----------|-------------|-------------|--------|----------|-----------|
|                            | **mean**      | **se_mean** | **sd** | **2.5%** | **97.5%** | **mean**    | **se_mean** | **sd** | **2.5%** | **97.5%** |
| &sigma;<sub>y</sub>        | 0.036         | 0.000       | 0.001  | 0.034    | 0.037     | 0.036       | 0.000       | 0.001  | 0.034    | 0.037     |
| &sigma;<sub>&mu;</sub>     | 0.038         | 0.000       | 0.001  | 0.035    | 0.040     | 0.038       | 0.000       | 0.001  | 0.035    | 0.040     |
| &phi;                      | 0.619         | 0.000       | 0.078  | 0.461    | 0.766     | 0.619       | 0.001       | 0.079  | 0.461    | 0.769     |
| *c*<sub>&alpha;</sub>      | 0.000         | 0.000       | 0.002  | -0.003   | 0.004     | 0.000       | 0.000       | 0.002  | -0.003   | 0.004     |
| &sigma;<sub>&alpha;</sub>  | 0.017         | 0.000       | 0.002  | 0.014    | 0.020     | 0.017       | 0.000       | 0.002  | 0.014    | 0.021     |
| &sigma;<sub>&lambda;</sub> | 0.463         | 0.002       | 0.080  | 0.316    | 0.631     | 0.463       | 0.002       | 0.081  | 0.317    | 0.632     |
| &sigma;<sub>&gamma</sub>   | 0.204         | 0.001       | 0.032  | 0.144    | 0.273     | 0.209       | 0.001       | 0.034  | 0.145    | 0.282     |

![image](https://user-images.githubusercontent.com/65618755/162422597-7db0f523-07d9-4d78-8d89-1c17360b0c91.png)

**Fig. 1.** Estimated latent processes &mu;<sub>i,t</sub> of the model with wide priors for some example sites.

![image](https://user-images.githubusercontent.com/65618755/162422740-54e7b75c-96ce-42c9-a36d-ae5eb7cb717a.png)

**Fig. 2.** Estimated trend &alpha;<sub>t</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/162422881-bf917a92-2102-482a-a91a-487d81e3bffa.png)

**Fig. 3.** Estimated trend intensity parameters &lambda;<sub>i</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/162423039-b6ea91b1-c8de-43ca-8142-68e583f907fd.png)

**Fig. 4.** Estimated short-term parameters &beta;<sub>i,j</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/162423235-a999f938-23ad-42b7-948b-f89a57c48899.png)

**Fig. 5.** Estimated error correction parameters &gamma;<sub>i,j,t</sub> of the model with wide priors.


### Comparison with classic SAR

A simpler SAR model was also fitted to the same data to compare the latent level outputs &mu;<sub>i,t</sub>. In this case the model did not include the trend &alpha;<sub>t</sub> or the error correction terms &gamma;<sub>i,j,t</sub>, and the short-term coefficients &beta;<sub>i,j</sub> were replaced with a single spatial dependency parameter, &rho;. More formally, we defined &mu;<sub>t</sub> - &mu;<sub>t-1</sub> = &rho;*W*(&mu;<sub>t</sub> - &mu;<sub>t-1</sub>) + &epsilon;<sub>t</sub>, leading to &mu;<sub>t</sub> &sim; N(&mu;<sub>t-1</sub>,&sigma;<sub>&mu;</sub><sup>2</sup>(*I* -&rho; *W*)<sup>-1</sup>((*I* - &rho;*W*)<sup>-1</sup>)<sup>T</sup>, where *W* is the fixed adjacency matrix after spectral normalization, and -1 < &rho; < 1 is the unknown spatial dependency parameter. Thus in addition to latent &mu;, this model contained only three unknown parameters, &sigma;<sub>y</sub>, &sigma;<sub>&mu;</sub>, and &rho;. The estimates of the &mu;<sub>i,t</sub> (plotted below) resembled those of our main model except that the posterior intervals were wider with reference to the main model.

**Table 2.** Scalar parameters of the SAR model.
|                      | mean  | se_mean | sd    | 2.5%  | 97.5% | n_eff    | Rhat  |
|----------------------|-------|---------|-------|-------|-------|----------|-------|
| &sigma;<sub>y</sub>  | 0.034 | 0.000   | 0.001 | 0.032 | 0.035 | 825.897  | 1.002 |
| &sigma;<sub>mu</sub> | 0.055 | 0.000   | 0.001 | 0.054 | 0.057 | 1185.511 | 1.002 |
| &rho;                | 0.778 | 0.000   | 0.012 | 0.755 | 0.801 | 3163.835 | 1.001 |


![image](https://user-images.githubusercontent.com/65618755/158409784-0a38869a-f857-4d67-a3ca-93ea2ed920f9.png)

**Fig. 6.** Estimated latent processes &mu;<sub>i,t</sub> of the classic SAR model for some example sites.

![image](https://user-images.githubusercontent.com/65618755/162441504-68d6822e-3981-42b4-b426-53789ba65f3f.png)

**Fig. 7.** Estimated latent processes &mu;<sub>i,t</sub> of the classic SAR model (blue line and ribbon) and the model in the paper (yellow line and ribbon) for some example sites.
