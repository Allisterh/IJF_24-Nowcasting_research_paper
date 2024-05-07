# IJF_24

Replication code for the paper Labonne, P. (2024). “Asymmetric
Uncertainty: Nowcasting using skewness in real-time data.”
*International Journal of Forecasting*

The score-driven models are written `C++` to make use of the automatic
differentiation library `CppAD`. This code is then wrapped into a
package using `Rcpp`. The package can be installed with:

``` r
devtools::install_github("paullabonne/IJF_24/RcppScoreDrivenDFM")
```

- `R/intro.R` loads/install all the packages needed.
- `R/sd_dfm.R` contains functions for estimation, forecast evaluation,
  simulation etc.
- The `quarto` files in the notebooks folder replicate the paper, i.e.,
  build the vintages, real-time recursive estimation, result analysis
  and simulations.
