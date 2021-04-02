
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GGMnonreg

[![CircleCI build
status](https://circleci.com/gh/donaldRwilliams/GGMnonreg.svg?style=svg)](https://circleci.com/gh/donaldRwilliams/GGMnonreg)

The goal of **GGMnonreg** is to estimate non-regularized Gaussian
graphical models. Note that the title is a bit of a misnomer, in that
Ising and mixed graphical models are also supported.

Graphical modeling is quite common in fields with *wide* data, that is,
when there are more variables than observations. Accordingly, many
regularization-based approaches have been developed for those kinds of
data.

More recently, graphical modeling has emerged in psychology (Epskamp et
al. 2018), where the data is typically long or low-dimensional \[*p* \<
*n*, Williams et al. (2019); williams\_rethinking\]. The primary purpose
of **GGMnonreg** is to provide methods specifically for low-dimensional
data (e.g., those common to psychopathology networks).

## Supported Models

  - Gaussian graphical model. The following data types are supported.
      - Gaussian
      - Ordinal
      - Binary
  - Ising model (Marsman et al. 2017)
  - Mixed graphical model

## Additional methods

The following are also included

  - Expected network replicability (Williams 2020)
  - Compare Gaussian graphical models
  - Measure of parameter uncertainty (Williams et al. 2019)
  - Edge inclusion “probabilities”
  - Network visualization with **ggplot2** (Wickham 2016)
  - Constrained precision matrix (the network, given an assumed graph)
  - Predictability (variance explained)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMnonreg")
```

## References

<div id="refs" class="references">

<div id="ref-Epskamp2018ggm">

Epskamp, Sacha, Lourens J. Waldorp, Rene Mottus, and Denny Borsboom.
2018. “The Gaussian Graphical Model in Cross-Sectional and Time-Series
Data.” *Multivariate Behavioral Research* 53 (4): 453–80.
<https://doi.org/10.1080/00273171.2018.1454823>.

</div>

<div id="ref-marsman_2018">

Marsman, M, D Borsboom, J Kruis, S Epskamp, R van Bork, L J Waldorp, By
Taylor, L Waldorp, H L J van der Maas, and G Maris. 2017. “An
Introduction to Network Psychometrics: Relating Ising Network Models to
Item Response Theory Models.” *Taylor & Francis* 53 (1): 15–35.
<https://doi.org/10.1080/00273171.2017.1379379>.

</div>

<div id="ref-ggplotpackage">

Wickham, Hadley. 2016. *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York. <http://ggplot2.org>.

</div>

<div id="ref-williams2020learning">

Williams, Donald R. 2020. “Learning to Live with Sampling Variability:
Expected Replicability in Partial Correlation Networks.” *PsyArXiv*.
<https://doi.org/10.31234/osf.io/fb4sa>.

</div>

<div id="ref-williams2019nonregularized">

Williams, Donald R, Mijke Rhemtulla, Anna C Wysocki, and Philippe Rast.
2019. “On Nonregularized Estimation of Psychological Networks.”
*Multivariate Behavioral Research* 54 (5): 719–50.
<https://doi.org/10.1080/00273171.2019.1575716>.

</div>

</div>
