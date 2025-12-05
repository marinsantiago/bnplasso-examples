# bnplasso-examples

This repository contains source code to reproduce the results from the article "[Adaptive Shrinkage with a Nonparametric Bayesian Lasso](https://doi.org/10.1080/10618600.2025.2572327)" (Marin et al., 2025+).

The `bnplasso` R package is available at the Github repository: [https://github.com/marinsantiago/bnplasso](https://github.com/marinsantiago/bnplasso).

If you wish to reproduce the results from Marin et al. (2025+), you should download the version of the 
package employed at that time (`bnplasso` v0.1.0). That can easily be done by running in R

``` r
# install.packages("devtools")
devtools::install_github("marinsantiago/bnplasso@3c87169")
```

Alternatively, if wish to reproduce the results from Marin et al. (2025+), you can also install
the package from the `bnplasso` folder in the supplementary materials to Marin et al. (2025+):

  1. In R, set your working directory to the folder `bnplasso`.
  
  2. Run the following R code:
  
``` r
# install.packages("devtools")
devtools::build()
devtools::install()
```

Detailed instructions on how to reproduce the results from the article are presented in the README file, located in the Supplementary Materials accompanying the main manuscript, available at [https://doi.org/10.1080/10618600.2025.2572327](https://doi.org/10.1080/10618600.2025.2572327).

## <a name="cite"></a> Citation

If you use any part of this code in your work, please consider citing our *JCGS* paper:

```TeX
@article{marin_bnplasso,
  title   = {Adaptive Shrinkage with a Nonparametric Bayesian Lasso},
  author  = {Santiago Marin and Bronwyn Loong and Anton H. Westveld},
  journal = {Journal of Computational and Graphical Statistics},
  year    = {2025},
  doi     = {10.1080/10618600.2025.2572327},
}
```

## <a name="refs"></a> References

Marin, S., Loong, B., and Westveld, A. H. (2025+), "Adaptive Shrinkage with a Nonparametric Bayesian Lasso." *Journal of Computational and Graphical Statistics*. [doi:10.1080/10618600.2025.2572327](https://doi.org/10.1080/10618600.2025.2572327)

</br>
