# spopt <a href="https://walker-data.com/spopt/"><img src="man/figures/logo.png" align="right" height="120" alt="spopt website" /></a>

The **spopt** R package provides R-native implementations of spatial optimization algorithms for regionalization, facility location, and market analysis. Inspired by [Python's PySAL spopt](https://pysal.org/spopt/), the package brings these powerful algorithms to R users with an sf-first API and a Rust backend for performance.

Install from r-universe (recommended):

```r
install.packages("spopt", repos = "https://walkerke.r-universe.dev")
```

Or install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("walkerke/spopt-r")
```

Read through these vignettes to learn how to use the package:

- [Getting started with __spopt__](https://walker-data.com/spopt/articles/getting-started.html)

- [Regionalization](https://walker-data.com/spopt/articles/regionalization.html)

- [Facility location](https://walker-data.com/spopt/articles/facility-location.html)

- [Huff model](https://walker-data.com/spopt/articles/huff-model.html)

- [Travel-time cost matrices](https://walker-data.com/spopt/articles/travel-time-matrices.html)

## Support and how to learn more

If you find this project useful in your work and would like to ensure continued development of the package, you can provide support in the following ways:

* [Chip in some funds to support package development via PayPal](https://www.paypal.com/paypalme/walkerdata/);
* Set up a consulting engagement or workshop through Walker Data to help you implement __spopt__ in your project. Send a note to <kyle@walker-data.com> if you are interested;
* File an issue - or even better, a pull request - at https://github.com/walkerke/spopt-r/issues.

To stay on top of package updates / new features and to get information about trainings, [be sure to sign up for the Walker Data mailing list here](https://walker-data.us15.list-manage.com/subscribe?u=1829a68a5eda3d301119fdcd6&id=c4a53d2961).
