---
title: 'GGMnonreg: Non-Regularized Gaussian Graphical Models'
tags:
- Graphical models
- partial correlations
- Mixed graphical model
- Ising model
date: "03 April 2021"
output: pdf_document
affiliations:
- name: Department of Psychology, University of California, Davis
index: 1
citation_author: Williams
authors:
- name: Donald R. Williams
year: 2021
bibliography: inst/REFERENCES.bib
affiliation: 1
---
  
# Summary
Studying complex relations in multivariate datasets is a common task across the sciences. Cognitive neuroscientists model brain connectivity with the goal of unearthing functional and structural associations between cortical regions [@ortiz_2015]. In clinical psychology, researchers wish to better understand the intricate web of symptom interrelations that underlie mental health disorders 
[@mcnally_2016; @borsboom_small_world]. To this end, graphical modeling has emerged as an oft-used tool in the chest of scientific inquiry. The basic idea is to characterize multivariate relations by learning the conditional dependence structure. The cortical regions or symptoms are *nodes* and the featured connections linking nodes are *edges* that graphically represent the conditional dependence structure.

Graphical modeling is quite common in fields with wide data, that is, when there are more variables ($p$) than observations ($n$). Accordingly, many regularization-based approaches have been developed for those kinds of data. There are key drawbacks of regularization, including, but not limited to, the fact that obtaining a valid measure of parameter uncertainty is very (very) difficult [@Buhlmann2014] and there can
be an inflated false positive rate [see for example, @williams2019nonregularized].

More recently, graphical modeling has emerged in psychology (Epskamp et al. 2018), where the data is typically long or low-dimensional ($p < n$; @williams2019nonregularized, @williams_rethinking). The primary purpose of GGMnonreg is to provide methods that were specifically designed for low-dimensional data (e.g., those common in the social-behavioral sciences).


# Statement of Need
The following were designed specifically for low-dimensional data, for which there is a dearth of methodology.

## Supported Models
* Gaussian graphical model. The following data types are supported.
  + Gaussian 
  + Ordinal
  + Binary
* Ising model [@marsman_2018]
* Mixed graphical model

## Additional methods
The following are also included

* Expected network replicability [@williams2020learning]
* Compare Gaussian graphical models
* Measure of parameter uncertainty [@williams2019nonregularized]
* Edge inclusion "probabilities" [e.g., Figure 6.4 in  @Hastie2015]
* Network visualization
* Constrained precision matrix [the network, given an assumed graph, see p. 631 in @hastie2009elements]
* Predictability [variance explained for each node, @Haslbeck2018]

## Implementation
```{r}
# make binary
Y <- ifelse(ptsd[,1:5] == 0, 0, 1)

# fit model
fit <- ising_search(Y, IC = "BIC", 
                    progress = FALSE)

fit
#>          1        2        3        4        5
#> 1 0.000000 1.439583 0.000000 1.273379 0.000000
#> 2 1.439583 0.000000 1.616511 0.000000 1.182281
#> 3 0.000000 1.616511 0.000000 1.716747 1.077322
#> 4 1.273379 0.000000 1.716747 0.000000 1.662550
#> 5 0.000000 1.182281 1.077322 1.662550 0.000000
```

Note the same code, more or less, is also used for GGMs and mixed graphical models.

## Network Visualization
A key aspect of graphical modeling is visualizing the conditional dependence structure. To this end, 
**GGMnonreg** makes network plots with **ggplot2** [@ggplotpackage].

```
plot(get_graph(fit), 
     node_names = colnames(Y), 
     edge_magnify = 2)
```
![Conditional Dependence Structure](man/figures/figure_1.png)

# Acknowledgements
DRW was supported by a National Science Foundation Graduate Research Fellowship
under Grant No. 1650042


# References

