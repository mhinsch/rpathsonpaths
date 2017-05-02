
Installing the package
======================

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

-   **`spread_dirichlet`**: simulate spread of genetic material on the network using a Dirichlet distribution to approximate genetic drift.

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

inp <- c(0, 0, 1, 2, 3, 1, 5)
outp <- c(1, 2, 3, 3, 4, 4, 2)
rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
edgelist <- data.frame(inputs=inp, outputs=outp, rates=rates)
ext <- data.frame(nodes=c(0, 5), rates=c(0.5, 0.5))
```

Now we can create a *popsnetwork* object:

``` r

netraw <- popsnetwork(edgelist, ext, 0.1)
netraw
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1
#> 1    0.55    1
#> 2    1.1 2
#> 3    0.357   0.6
#> 4    0.695   1.1
#> 5    0.5 1
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.5
#> 0    2   1.5 0.75
#> 1    3   0.5 0.275
#> 2    3   0.1 0.055
#> 3    4   1   0.595
#> 1    4   0.1 0.055
#> 5    2   0.5 0.25
plot(netraw)
#> Loading required namespace: igraph
```

![](figs/popsnetwork-1.png)

Injecting initial populations and simulating spread
---------------------------------------------------

We set up the initial populations of the source nodes:

``` r
freqs <- matrix(c(0.1, 0.4, 0.3, 0.1, 0.2, 0.2, 0.4, 0.3), ncol=4, nrow=2)
netini <- set_allele_freqs(netraw, list(nodes=c(0, 5), frequencies=freqs))
netini
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1   0.1 0.3 0.2 0.4
#> 1    0.55    1
#> 2    1.1 2
#> 3    0.357   0.6
#> 4    0.695   1.1
#> 5    0.5 1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.5
#> 0    2   1.5 0.75
#> 1    3   0.5 0.275
#> 2    3   0.1 0.055
#> 3    4   1   0.595
#> 1    4   0.1 0.055
#> 5    2   0.5 0.25
```

Now we are ready to simulate spread of genetic material in the network. Note that every call to spread\_dirichlet will generally produce a different allele frequency distribution.

``` r
netdir1 <- spread_dirichlet(netini, 1.7)
netdir1
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1   0.1 0.3 0.2 0.4
#> 1    0.55    1   0.336698    0.219971    0.0142882   0.429043
#> 2    1.1 2   0.226885    0.54888 0.120218    0.104017
#> 3    0.357   0.6 0.898183    0.0764407   6.26606e-05 0.0253133
#> 4    0.695   1.1 0.834274    0.0898766   0.0548936   0.0209557
#> 5    0.5 1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.5
#> 0    2   1.5 0.75
#> 1    3   0.5 0.275
#> 2    3   0.1 0.055
#> 3    4   1   0.595
#> 1    4   0.1 0.055
#> 5    2   0.5 0.25
netdir2 <- spread_dirichlet(netini, 1.7)
netdir2
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1   0.1 0.3 0.2 0.4
#> 1    0.55    1   0.00180906  0.0759773   0.438173    0.484041
#> 2    1.1 2   0.105474    0.159487    0.529805    0.205233
#> 3    0.357   0.6 0.0355106   0.00225909  0.195762    0.766468
#> 4    0.695   1.1 1.21519e-05 3.50116e-06 0.110928    0.889057
#> 5    0.5 1   0.4 0.1 0.2 0.3
#> 
#> Links:
#> 
#> from to  rate    infected
#> 0    1   1   0.5
#> 0    2   1.5 0.75
#> 1    3   0.5 0.275
#> 2    3   0.1 0.055
#> 3    4   1   0.595
#> 1    4   0.1 0.055
#> 5    2   0.5 0.25
```

Drawing samples
---------------

Finally we can draw samples from our simulated population.

``` r
samplconf <- data.frame(nodes=c(2, 4), N=c(10, 10))
draw_isolates(netdir1, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        1        6        2        1
#> 2    4        9        0        1        0
draw_isolates(netdir2, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        2        2        5        1
#> 2    4        0        0        1        9
```

Contributors (in alphabetic order):
===================================

-   [Martin Hinsch](https://github.com/mhinsch)
-   [Thibaut Jombart](https://github.com/thibautjombart)

See details of contributions on: <br> <https://github.com/mhinsch/rpathsonpaths/graphs/contributors>

Contributions are welcome via **pull requests**.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

**Maintainer:** Martin Hinsch (<hinsch.martin@gmail.com>)
