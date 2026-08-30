
# LUCIDus

<!-- badges: start -->

<!-- badges: end -->

LUCIDus implements **LUCID** (Latent Unknown Clusters by Integrating
multi-omics Data): a model that finds latent subgroups of subjects that
are simultaneously predictable from a set of exposures, characterized by
distinct multi-omics profiles, and associated with a health outcome. It
supports three integration strategies – **early**, **parallel**, and
**serial** – along with feature selection, missing-data handling,
prediction, g-computation, and bootstrap inference.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/LUCIDus")
```

## A minimal example: parallel integration on the bundled HELIX data

This walks through fitting a LUCID **parallel** model – one latent
cluster variable per omics layer, fit jointly – on
`simulated_HELIX_data`, the simulated HELIX-study dataset bundled with
the package. Parallel integration is the right choice here because the
omics data comes as three genuinely distinct layers (methylome,
transcriptome, miRNA); rather than pooling them into one matrix, each
layer gets its own cluster variable, estimated together. For a much more
detailed, narrated walkthrough of all three architectures (including the
two-step screen-then-refit workflow for real feature selection), see
`vignette("lucid_3models_normal_outcome")`.

``` r
library(LUCIDus)

# The bundled HELIX dataset: one exposure (postnatal mercury exposure),
# three omics layers, two covariates, and one outcome (a liver enzyme).
data(simulated_HELIX_data)
d <- simulated_HELIX_data

G <- as.matrix(d$phenotype$hs_hg_m_scaled)
colnames(G) <- "mercury_exposure"

Z <- list(
  methylome     = d$methylome,
  transcriptome = d$transcriptome,
  miRNA         = d$miRNA
)

Co <- data.frame(
  age = d$phenotype$hs_child_age_yrs_None,
  sex = as.numeric(d$phenotype$e3_sex_None == "male")
)

Y <- as.matrix(d$phenotype$ck18_scaled)
colnames(Y) <- "ck18"

fit <- estimate_lucid(
  G = G, Z = Z, Y = Y, CoG = Co, CoY = Co,
  lucid_model = "parallel",
  family = "normal",
  K = c(2, 2, 2),
  seed = 2026
)
#> Fitting LUCID parallel model (3 layers)...
#> Finished LUCID parallel model.

class(fit)
#> [1] "lucid_parallel"
```

`lucid_model = "parallel"` selects the architecture; `K = c(2, 2, 2)`
asks for 2 clusters in each of the three layers (in `Z`’s order).
`CoG`/`CoY` are covariates adjusted for in the exposure-to-cluster and
cluster-to-outcome models respectively – here, the same age/sex
covariates for both.

### Reading the fitted model

``` r
summary(fit)
#> 
#> ====================================================
#> LUCID Parallel: Model Summary
#> ====================================================
#> 
#> Model specification
#>   Family                 : gaussian
#>   Number of observations : 420
#>   Clusters per layer     : 2, 2, 2
#> 
#> Missing-data profile by layer
#>   Layer 1 listwise rows : 0 / 420 (0.0%)
#>   Layer 1 sporadic rows : 0 / 420 (0.0%)
#>   Layer 1 missing cells : 0 / 4200 (0.0%)
#>   Layer 2 listwise rows : 0 / 420 (0.0%)
#>   Layer 2 sporadic rows : 0 / 420 (0.0%)
#>   Layer 2 missing cells : 0 / 4200 (0.0%)
#>   Layer 3 listwise rows : 0 / 420 (0.0%)
#>   Layer 3 sporadic rows : 0 / 420 (0.0%)
#>   Layer 3 missing cells : 0 / 4200 (0.0%)
#> 
#> Feature selection overview
#>   G features selected    : 1 / 1 (100.0%)
#>   G features by layer
#>     Layer 1              : 1 / 1 (100.0%)
#>     Layer 2              : 1 / 1 (100.0%)
#>     Layer 3              : 1 / 1 (100.0%)
#>   Z features
#>     Layer 1 selected     : 10 / 10 (100.0%)
#>     Layer 1 multi-cluster: 10
#>     Layer 2 selected     : 10 / 10 (100.0%)
#>     Layer 2 multi-cluster: 10
#>     Layer 3 selected     : 10 / 10 (100.0%)
#>     Layer 3 multi-cluster: 10
#> 
#> Model fit statistics
#>   Log-likelihood         : -16169.04
#>   BIC                    : 34808.55
#>   Number of parameters   : 409
#> 
#> Regularization
#>   Rho_G                  : 0.000
#>   Rho_Z_Mu               : 0.000
#>   Rho_Z_Cov              : 0.000
#> 
#> Detailed parameter estimates
#> (1) Y (continuous outcome): intercept, effects of each non-reference latent cluster for each layer of Y (and effect of covariates if included) 
#>                   Gamma
#> (Intercept) -0.82350744
#> Layer1_LC2   0.95020537
#> Layer2_LC2   0.19542997
#> Layer3_LC2   0.47724613
#> age         -0.00863345
#> sex          0.02341732
#> 
#> (2) Z: mean of omics data for each latent cluster of each layer 
#> Layer 1
#> 
#>               mu_cluster1 mu_cluster2
#> cg_GRHL3       0.18408265 -0.36284441
#> cg_BTF3L4      0.06723869 -0.17541351
#> cg_AL358472.7  0.17008486 -0.36863454
#> cg_HDGF        0.21746220 -0.36714732
#> cg_TDRD5       0.25376985 -0.36018991
#> cg_CSRNP3      0.08050203 -0.09532648
#> cg_HSPD1       0.17490777 -0.26096784
#> cg_EPM2AIP1   -0.13860873  0.34406638
#> cg_AC025171.1  0.34470233 -0.44257881
#> cg_VTRNA1_3   -0.11214165  0.08113805
#> 
#> Layer 2
#> 
#>                   mu_cluster1 mu_cluster2
#> tc_TC01006069_nc -0.522270317  0.47930061
#> tc_SLC9A4        -0.006309301  0.01450134
#> tc_RAB6C_AS1     -0.041387945 -0.02026131
#> tc_LOC100129029   0.153725322 -0.04718847
#> tc_BRE           -0.088592471  0.10032251
#> tc_TC03001220_nc -0.361877314  0.28246497
#> tc_TC04002114_nc -0.128442683  0.11614619
#> tc_TC04002369_nc -0.375653735  0.52775346
#> tc_BEND4         -0.255524235  0.16054010
#> tc_SLC9A3         0.202163670 -0.03271405
#> 
#> Layer 3
#> 
#>                mu_cluster1  mu_cluster2
#> miR.101.3p     0.031111557  0.064610758
#> miR.125a.5p   -0.153728401  0.079688810
#> miR.125b.1.3p  0.039869641  0.040503869
#> miR.127.3p    -0.415405400  0.314634153
#> miR.140.5p     0.312852637 -0.009306536
#> miR.142.3p    -0.102400018  0.086767597
#> miR.144.5p    -0.005121102  0.027388439
#> miR.19a.3p    -0.276896191  0.165091012
#> miR.19b.3p    -0.126646333  0.128221532
#> miR.21.5p      0.003262879  0.055314865
#> 
#> (3) E: intercept and odds ratio of being assigned to each latent cluster for each exposure for each layer 
#> Layer 1
#> 
#>                                beta        OR
#> (Intercept).cluster2      2.4308007 11.367981
#> mercury_exposure.cluster2 0.9694861  2.636589
#> 
#> Layer 2
#> 
#>                                 beta       OR
#> (Intercept).cluster2      -1.0848810 0.337942
#> mercury_exposure.cluster2  0.4478715 1.564978
#> 
#> Layer 3
#> 
#>                                 beta        OR
#> (Intercept).cluster2      -0.1489437 0.8616176
#> mercury_exposure.cluster2 -0.0522552 0.9490866
```

`summary()` reports, per layer: the model specification, the exposure’s
effect on cluster membership (as a log-odds/odds ratio), and each
cluster’s mean omics profile. Two extractor functions read the same
information programmatically, auto-detecting that `fit` is a parallel
model with no further input:

``` r
# Each subject's hard cluster assignment, one vector per layer.
str(get_cluster_assignment(fit))
#> List of 3
#>  $ : num [1:420] 1 1 1 2 1 1 1 1 2 2 ...
#>  $ : num [1:420] 2 2 2 1 2 1 2 2 2 2 ...
#>  $ : num [1:420] 2 1 2 1 2 1 1 1 2 2 ...

# The 5 omics features that most separate the clusters, one layer's worth.
get_top_omics_features(fit, top_n = 5)[["methylome"]]
#> cg_AC025171.1      cg_TDRD5       cg_HDGF cg_AL358472.7      cg_GRHL3 
#>     0.5750393     0.4513173     0.4322713     0.3922378     0.3803111
```

### Visualizing the clusters’ omics profiles

`plot_cluster_omic_profile()` shows which omics features distinguish the
clusters, and in which direction – one figure per layer for a parallel
fit:

``` r
plot_cluster_omic_profile(fit, top_n = 8)[["methylome"]]
```

![](man/figures/README-profile-1.png)<!-- -->

## Learn more

Once installed, every vignette below is available with
`vignette("<name>")`. Rendered copies are also in [`docs/`](docs/) in
this repository, viewable directly in a browser without installing the
package:

- [`lucid_3models_normal_outcome`](https://htmlpreview.github.io/?https://github.com/ContiLab-usc/LUCIDus-3.0/blob/main/docs/lucid_3models_normal_outcome.html)
  – a full early/parallel/serial walkthrough, including feature
  selection and bootstrap inference, for a **continuous** outcome.
- [`lucid_3models_binary_outcome`](https://htmlpreview.github.io/?https://github.com/ContiLab-usc/LUCIDus-3.0/blob/main/docs/lucid_3models_binary_outcome.html)
  – the same walkthrough for a **binary** outcome.
- [`lucidus_full_functionality_guide`](https://htmlpreview.github.io/?https://github.com/ContiLab-usc/LUCIDus-3.0/blob/main/docs/lucidus_full_functionality_guide.html)
  – a breadth-first tour of every exported function.
- `NEWS.md` – what changed in the current release.
