
[![Travis-CI Build Status](https://travis-ci.org/mhinsch/rpathsonpaths.svg?branch=master)](https://travis-ci.org/mhinsch/rpathsonpaths)

Installing the package
======================

Prerequisites:

-   C++ boost libraries

-   *Rcpp* (&gt;= 0.12.9)

Some of the functionality will also require the *igraph* package.

*rpathsonpaths* is not on CRAN yet, therefore for now please install the development version of the package from *github*:

``` r
devtools::install_github("mhinsch/rpathsonpaths")
```

Note that this requires the package *devtools* installed.

What does it do?
================

*rpathsonpaths* provides functions to simulate the spread of a population in a substrate that itself is transported across a directed graph. The paradigmatic example is the spread of pathogens in food transport networks, but other instances might include for examples infections in migrating animals.

The most important functions include:

-   **`popsnetwork`**: create a network object from an edgelist (describing the transport graph) and a list of external sources.

-   **`set_allele_freq`**: set the population composition in source nodes.

-   **`popgen_dirichlet`**: simulate spread of genetic material on the network using a Dirichlet distribution to approximate genetic drift.

-   **`draw_isolates`**: draw a number of random individuals from a set of nodes using the simulated allele frequencies.

-   **`plot`**: plot a *popsnetwork* object (currently uses igraph, but more backends are planned).

Resources
=========

Vignettes
---------

An overview of *rpathsonpaths* is provided below in the worked example below. More detailed tutorials are distributed as vignettes with the package:

``` r
vignette("overview", package="rpathsonpaths")
```

Websites
--------

The following websites are available:

-   The *rpathsonpaths* project on *github*, useful for developers, contributors, and users wanting to post issues, bug reports and feature requests: <br> [http://github.com/mhinsch/rpathsonpaths](http://github.com/reconhub/rpathsonpaths)

Getting help online
-------------------

Bug reports and feature requests should be posted on *github* using the [*issue*](http://github.com/mhinsch/rpathsonpaths/issues) system.

A quick overview
================

The following worked example provides a brief overview of the package's functionalities. See the [*vignettes section*](#vignettes) for more detailed tutorials.

Preparing the network
---------------------

First, we prepare an edge list describing the network, and the initial state of the source nodes:

``` r
library(rpathsonpaths)

inp <- c(0L, 0L, 1L, 2L, 3L, 1L, 5L)
outp <- c(1L, 2L, 3L, 3L, 4L, 4L, 2L)
rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
edgelist <- data.frame(inp, outp, rates)
ext <- data.frame(c(0L, 5L), c(0.5, 0.5))
```

Now we can create a *popsnetwork* object:

``` r

netraw <- popsnetwork(edgelist, ext, 0.1)
netraw
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.55    1
#> 1    0.595   1
#> 2    0.595   2
#> 3    0.6355  0.6
#> 4    0.668636    1.1
#> 5    0.55    1
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.55
#> 0    2   1.5 0.55
#> 1    3   0.5 0.595
#> 2    3   0.1 0.595
#> 3    4   1   0.6355
#> 1    4   0.1 0.595
#> 5    2   0.5 0.55
plot(netraw)
#> Loading required namespace: igraph
```

![](figs/popsnetwork-1.png)

Injecting initial populations and simulating spread
---------------------------------------------------

We set up the initial populations of the source nodes:

``` r
freqs <- matrix(c(0.1, 0.4, 0.3, 0.1, 0.2, 0.2, 0.4, 0.3), ncol=4, nrow=2)
netini <- set_allele_freqs(netraw, list(c(0L, 5L), freqs))
netini
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.55    1   0.1 0.3 0.2 0.4
#> 1    0.595   1   0   0   0   0
#> 2    0.595   2   0   0   0   0
#> 3    0.6355  0.6 0   0   0   0
#> 4    0.668636    1.1 0   0   0   0
#> 5    0.55    1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.55
#> 0    2   1.5 0.55
#> 1    3   0.5 0.595
#> 2    3   0.1 0.595
#> 3    4   1   0.6355
#> 1    4   0.1 0.595
#> 5    2   0.5 0.55
```

Now we are ready to simulate spread of genetic material in the network. Note that every call to popgen\_dirichlet will generally produce a different allele frequency distribution.

``` r
netdir1 <- popgen_dirichlet(netini, 1.7)
netdir1
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.55    1   0.1 0.3 0.2 0.4
#> 1    0.595   1   0.00162633  0.0995237   0.508113    0.390737
#> 2    0.595   2   0.25112 0.620695    0.00135567  0.12683
#> 3    0.6355  0.6 0.0758687   0.0889155   0.0142789   0.820937
#> 4    0.668636    1.1 0.0010684   0.25515 0.0142712   0.72951
#> 5    0.55    1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.55
#> 0    2   1.5 0.55
#> 1    3   0.5 0.595
#> 2    3   0.1 0.595
#> 3    4   1   0.6355
#> 1    4   0.1 0.595
#> 5    2   0.5 0.55
netdir2 <- popgen_dirichlet(netini, 1.7)
netdir2
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.55    1   0.1 0.3 0.2 0.4
#> 1    0.595   1   0.0257071   0.00527239  0.00192251  0.967098
#> 2    0.595   2   0.157011    0.51551 0.0656875   0.261791
#> 3    0.6355  0.6 0.0226906   0.104864    6.00565e-11 0.872445
#> 4    0.668636    1.1 1.83311e-12 0.0123462   7.32917e-07 0.987653
#> 5    0.55    1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.55
#> 0    2   1.5 0.55
#> 1    3   0.5 0.595
#> 2    3   0.1 0.595
#> 3    4   1   0.6355
#> 1    4   0.1 0.595
#> 5    2   0.5 0.55
```

Drawing samples
---------------

Finally we can draw samples from our simulated population.

``` r
samplconf <- data.frame(nodes=c(2L, 4L), N=c(10L, 10L))
draw_isolates(netdir1, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        1        8        0        1
#> 2    4        0        3        0        7
draw_isolates(netdir2, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        0        8        0        2
#> 2    4        0        0        0       10
```

Contributors (in alphabetic order):
===================================

-   [Martin Hinsch](https://github.com/mhinsch)
-   [Thibaut Jombart](https://github.com/thibautjombart)

See details of contributions on: <br> <https://github.com/mhinsch/rpathsonpaths/graphs/contributors>

Contributions are welcome via **pull requests**.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

**Maintainer:** Martin Hinsch (<hinsch.martin@gmail.com>)
