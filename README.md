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
| &sigma;<sub>y</sub>        | 0.036         | 0.000       | 0.001  | 0.034    | 0.038     | 0.036       | 0.000       | 0.001  | 0.034    | 0.038     |
| &sigma;<sub>mu</sub>       | 0.038         | 0.000       | 0.001  | 0.035    | 0.040     | 0.037       | 0.000       | 0.001  | 0.034    | 0.039     |
| &phi;                      | 0.620         | 0.001       | 0.078  | 0.462    | 0.769     | 0.622       | 0.001       | 0.077  | 0.469    | 0.770     |
| *c*<sub>a</sub>            | 0.000         | 0.000       | 0.002  | -0.003   | 0.004     | 0.000       | 0.000       | 0.002  | -0.003   | 0.004     |
| &sigma;<sub>a</sub>        | 0.017         | 0.000       | 0.002  | 0.014    | 0.021     | 0.016       | 0.000       | 0.002  | 0.013    | 0.020     |
| &sigma;<sub>&lambda;</sub> | 0.422         | 0.002       | 0.067  | 0.302    | 0.565     | 0.456       | 0.002       | 0.071  | 0.325    | 0.604     |
| &sigma;<sub>c</sub>        | 0.203         | 0.001       | 0.032  | 0.144    | 0.273     | 0.213       | 0.001       | 0.032  | 0.152    | 0.279     |

![image](https://user-images.githubusercontent.com/65618755/158409579-54d4b0fa-12c8-4577-b883-ab6f9fcf83d3.png)

**Fig. 1.** Estimated latent processes &mu;<sub>i,t</sub> of the model with wide priors for some example sites.

![image](https://user-images.githubusercontent.com/65618755/158410914-671cb854-6314-401a-ac7b-accd0dd0bb51.png)

**Fig. 2.** Estimated trend &alpha;<sub>t</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/158410226-681157bc-4a5b-430b-a858-5a59528974aa.png)

**Fig. 3.** Estimated trend intensity parameters &lambda;<sub>i</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/158410515-af2d9ead-fd24-4708-99f2-478f320189d0.png)

**Fig. 4.** Estimated short-term parameters &beta;<sub>i,j</sub> of the model with wide priors.

![image](https://user-images.githubusercontent.com/65618755/158410596-27642202-2f28-43e0-b7dc-e842cfb1f057.png)

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
