Bayesian Estimation of Hierarchical Linear Models from Incomplete Data:
Cluster-Level Interaction Effects
================
Dongho Shin

## GSExact Package

This repository provides the code for GSExact pacakge to perform
Bayesian estimation of Hierarchical Model (HLM) that is described in
paper “Bayesian Estimation of Hierarchical Linear Models from Incomplete
Data: Cluster-Level Interaction Effects and Small Sample Sizes” by Shin
et al (2024).

## Limitation

The current version of the package is limited to handling the
hierarchical linear model (HLM) presented in our manuscript:

$Y_{ij}=\beta_{0}+\beta_{1}C_{1j}+\beta_{2}C_{2j}+\beta_{3}X_{j}+\beta_{4}C_{1j}C_{2j}+u_{j}+e_{ij}$,

where $Y_{ij}$ is the outcome variable, $C_{1j}$ and $C_{2j}$ are
partially observed cluster-level variables with interaction effects, and
$X_{j}$ is a fully observed cluster-level variable. Additionally,
$u_j \sim \mathcal{N}(0, \tau)$ represents a level-2 unit-specific
random effect, and $e_{ij} \sim \mathcal{N}(0, \sigma^2)$ represents a
level-1 unit-specific random effect, both of which are independent.

Currently, the package assumes that the cluster size $n$ is consistent
across all clusters $J$.

This limitation primarily arises from the need to ensure compatibility.
For each unique HLM of interest, exact posterior distributions (for
parameters and missing covariates) must be derived to guarantee
compatible imputation models. These imputation models must then be
seamlessly integrated into the steps of the Gibbs sampler.

We are actively developing a more comprehensive software solution to
automate these complexities for users. Our aim is to create a versatile
tool capable of accommodating any desired HLM specification. This
ongoing development is expected to significantly enhance the
accessibility and applicability of our method.

In the meantime, researchers who wish to apply our method to HLMs
differing from the one presented in the manuscript are encouraged to
contact us at <shindh823@gmail.com>. We are happy to provide guidance on
adapting the current implementation to meet their specific needs.

### Dependencies

The GSExact package depends on the following packages: <br> lme4 <br>
MASS <br> monomvn <br> invgamma <br>

You can install these packages using:

``` r
install.packages("lme4")
install.packages("MASS")
install.packages("monomvn")
install.packages("invgamma")
```

### Installation

To install the GSExact package, use the following command:

``` r
install.packages("devtools")
library(devtools)
install_github("shind10/GSExact")
```

### Example

Here is an example of how to use the GSExact_basic function in the
GSExact package with the example data provided.

``` r
library(GSExact)
load("Example01.Rdata") #dataset name savedat
output <- GSExact_basic(input_data=example_data, J=36, n=4, burn=1000, iter=2500)
```

### Ouputs

The output consists of four lists: <br> output\[\[1\]\]: The estimated
parameters (calculated as the mean of post burn-in iterations). <br>
output\[\[2\]\]: The variance (ASE). <br> output\[\[3\]\]: The 95%
Bayesian credible interval.<br> output\[\[4\]\]: All estimated
parameters from every iteration (allowing you to compute the median or
generate a trace plot).
