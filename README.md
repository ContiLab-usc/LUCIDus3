
# LUCIDus Version 3: LUCID with Multiple Omics Data 
<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/LUCIDus?color=green)](https://cran.r-project.org/package=LUCIDus)
![](https://cranlogs.r-pkg.org/badges/grand-total/LUCIDus?color=blue)
[![](https://raw.githubusercontent.com/USCbiostats/badges/master/tommy-image-badge.svg)](https://image.usc.edu)
<!-- badges: end -->

This repository contains the source code for LUCIDus Version 3.0.3. The **LUCIDus** package implements the statistical method LUCID originally proposed in the research paper [A Latent Unknown Clustering Integrating
Multi-Omics Data (LUCID) with Phenotypic Traits](https://doi.org/10.1093/bioinformatics/btz667)
(*Bioinformatics*, 2020). The original LUCID conducts integrated clustering by using multi-view data, including genetic or environmental exposures, single-layer omics data, and  with/without outcome. **LUCIDus** features variable selection, incorporating missingness in omics data, visualization of the LUCID model via Sankey diagram, bootstrap inference, and functions for tuning model parameters. The missing data imputation mechanism was introduced in [An extension of latent unknown clustering integrating multi-omics data (LUCID) incorporating incomplete omics data](https://doi.org/10.1093/bioadv/vbae123)(*Bioinformatics Advances*,2024). We then extended the single-omic LUCID model to multi-omics in this research paper [LUCIDus: An R Package For Implementing Latent Unknown Clustering By Integrating Multi-omics Data (LUCID) With Phenotypic Traits](https://journal.r-project.org/articles/RJ-2024-012/)(*R Journal*,2024).

 **LUCID version 3**, a major update and enhancement from the original release, implements different integration strategies for multi-omics data with multiple layers, including LUCID early integration (early integration), LUCID in parallel (intermediate integration), and LUCID in serial (late integration). The following DAG illustrates the three different LUCID models for three integration strategies.

<p align="center">
  <img src="./figure/fig1-1.png" width="600">
</p>

If you are interested in the integration of omic data to estimate mediator or latent structures, please check out [Conti
Lab](https://contilab.usc.edu/about/) to learn more.



## Installation

You can install the development version of LUCIDus 3.0.3 from R CRAN with:

``` r
install.packages("LUCIDus")
```


## Workflow
The following figure illustrates the workflow of LUCIDus 3.0.3.
<p align="center">
  <img src="./figure/fig2-1.png" width="600">
</p>


## Usage

Please refer to Section 3 of the [R Journal article](https://journal.r-project.org/articles/RJ-2024-012/) for a full workflow illustration with real data.

Below are minimal examples showing the simplest usage of `lucid()` under **early**, **intermediate**, and **late** integrations.  
---

### Minimal Example: Early Integration 

```r
library(LUCIDus)

set.seed(123)

# Exposure (G): 1 variable
G <- matrix(rnorm(100), ncol = 1)

# Omics data (Z): 5 features
Z <- matrix(rnorm(100 * 5), ncol = 5)

# Outcome (Y): continuous
Y <- rnorm(100)

# Fit LUCID with early integration
fit_early <- lucid(
  G = G,
  Z = Z,
  Y = Y,
  family = "normal",
  K = 2
)

summary(fit_early)
```

### Minimal Example: Parallel Integration

---
```r
# Exposure
G <- matrix(rnorm(100), ncol = 1)

# Two omics layers (Z1 and Z2)
Z1 <- matrix(rnorm(100 * 4), ncol = 4)
Z2 <- matrix(rnorm(100 * 3), ncol = 3)

Z_list <- list(Z1 = Z1, Z2 = Z2)

# Outcome
Y <- rnorm(100)

fit_parallel <- lucid(
  G = G,
  Z = Z_list,
  Y = Y,
  family = "normal",
  lucid_model = "parallel",
  K = list(2,2)
)

summary(fit_parallel)
```

### Minimal Example: Serial Integration

---
```r
# Exposure
G <- matrix(rnorm(100), ncol = 1)

# Two omics layers for the serial pathway
Z1 <- matrix(rnorm(100 * 4), ncol = 4)
Z2 <- matrix(rnorm(100 * 3), ncol = 3)

Z_list <- list(Z1 = Z1, Z2 = Z2)

# Outcome
Y <- rnorm(100)

fit_serial <- lucid(
  G = G,
  Z = Z_list,
  Y = Y,
  family = "normal",
  lucid_model = "serial",
  K = list(2,2)
)

summary(fit_serial)
```

## Citation

If you use **LUCIDus**, please cite the following papers:

---

### **1. Original LUCID Method**

Peng C., Wang J., Asante I., Louie S., Jin R., Chatzi L., Casey G., Thomas D.C., Conti D.V. (2019).  
**A latent unknown clustering integrating multi-omics data (LUCID) with phenotypic traits.**  
*Bioinformatics*, btz667.  
https://doi.org/10.1093/bioinformatics/btz667

---

### **2. LUCIDus Package (R Journal Article for the Extension to Multi-omics)**

Zhao Y., Jia Q., Goodrich J.A., Conti D.V. (2024).  
**LUCIDus: An R Package for Implementing Latent Unknown Clustering by Integrating Multi-omics Data (LUCID) With Phenotypic Traits.**  
*R Journal*, 16(2), 2024.  
https://journal.r-project.org/articles/RJ-2024-012/

---

### **3. LUCID with Missing-Data Imputation (If the Missing Data Imputation Mechanism is Used)**

Zhao Y., Jia Q., Goodrich J., Darst B., Conti D.V. (2024).  
**An extension of latent unknown clustering integrating multi-omics data (LUCID) incorporating incomplete omics data.**  
*Bioinformatics Advances*, 4(1), vbae123.  
https://doi.org/10.1093/bioadv/vbae123

---

### **BibTeX**

You may obtain BibTeX entries in R using:

```r
print(citation("LUCIDus"), bibtex = TRUE)
toBibtex(citation("LUCIDus"))
options(citation.bibtex.max = 999)

