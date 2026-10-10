# spopt

The **spopt** R package provides R-native implementations of spatial
optimization algorithms for regionalization, facility location, route
optimization, and corridor routing. Inspired by [Python’s PySAL
spopt](https://pysal.org/spopt/), the package brings these algorithms to
R users with an sf-first API, a Rust backend for graph and routing
algorithms, and the [HiGHS](https://cran.r-project.org/package=highs)
solver for facility location.

Install the official release from CRAN:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"spopt"``)`

Or, install the development version from r-universe:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``'spopt'``, repos ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``'https://walkerke.r-universe.dev'``, ``'https://cloud.r-project.org'``)``)`

For developers who need to build from source (requires Rust), see
[CONTRIBUTING.md](https://github.com/walkerke/spopt-r/blob/main/CONTRIBUTING.md).

Read through these vignettes to learn how to use the package:

- [Getting started with
  **spopt**](https://walker-data.com/spopt-r/articles/getting-started.html)

- [Regionalization](https://walker-data.com/spopt-r/articles/regionalization.html)

- [Facility location
  problems](https://walker-data.com/spopt-r/articles/facility-location.html)

- [Huff model](https://walker-data.com/spopt-r/articles/huff-model.html)

- [Travel-time cost
  matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.html)

- [Routing
  optimization](https://walker-data.com/spopt-r/articles/routing.html)

- [Corridor routing on cost
  rasters](https://walker-data.com/spopt-r/articles/corridor-routing.html)

## Support and how to learn more

If you find this project useful in your work and would like to ensure
continued development of the package, you can provide support in the
following ways:

- [Chip in some funds to support package development via
  PayPal](https://www.paypal.com/paypalme/walkerdata/);
- Set up a consulting engagement or workshop through Walker Data to help
  you implement **spopt** in your project. Send a note to
  <kyle@walker-data.com> if you are interested;
- File an issue - or even better, a pull request - at
  <https://github.com/walkerke/spopt-r/issues>.

To stay on top of package updates / new features and to get information
about trainings, [be sure to sign up for the Walker Data mailing list
here](https://walker-data.us15.list-manage.com/subscribe?u=1829a68a5eda3d301119fdcd6&id=c4a53d2961).
