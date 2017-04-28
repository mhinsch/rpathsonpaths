
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

-   **`PopsNetwork`**: create a network object from an edgelist (describing the transport graph) and a list of external sources.

-   **`setAlleleFreq`**: set the population composition in source nodes.

-   **`spreadDirichlet`**: simulate spread of genetic material on the network using a Dirichlet distribution to approximate genetic drift.

-   **`drawIsolates`**: draw a number of random individuals from a set of nodes using the simulated allele frequencies.

-   **`plot`**: plot a *PopsNetwork* object (currently uses igraph, but more backends are planned).

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

Now we can create a *PopsNetwork* object:

``` r

netraw <- PopsNetwork(edgelist, ext, 0.1)
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
netini <- setAlleleFreqs(netraw, list(nodes=c(0, 5), frequencies=freqs))
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

Now we are ready to simulate spread of genetic material in the network. Note that every call to spreadDirichlet will generally produce a different allele frequency distribution.

``` r
netdir1 <- spreadDirichlet(netini, 1.7)
netdir1
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1   0.1 0.3 0.2 0.4
#> 1    0.55    1   0.0475387   1.33323e-05 0.111265    0.841183
#> 2    1.1 2   0.242394    0.0529073   0.530861    0.173838
#> 3    0.357   0.6 0.0148594   4.43345e-05 0.159789    0.825307
#> 4    0.695   1.1 6.71278e-20 0   0.0168881   0.983112
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
netdir2 <- spreadDirichlet(netini, 1.7)
netdir2
#> Nodes:
#> 
#> id   infected    input   alleles...
#> 0    0.5 1   0.1 0.3 0.2 0.4
#> 1    0.55    1   1.94912e-10 0.666281    0.198581    0.135137
#> 2    1.1 2   0.358183    0.125056    0.232029    0.284732
#> 3    0.357   0.6 0.0722427   0.432794    0.415673    0.0792901
#> 4    0.695   1.1 0.563539    0.0656922   0.0624262   0.308343
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
drawIsolates(netdir1, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        3        0        5        2
#> 2    4        0        0        1        9
drawIsolates(netdir2, samplconf)
#>   node allele_0 allele_1 allele_2 allele_3
#> 1    2        3        1        4        2
#> 2    4        6        0        0        4
```

Contributors (in alphabetic order):
===================================

-   [Martin Hinsch](https://github.com/mhinsch)
-   [Thibaut Jombart](https://github.com/thibautjombart)

See details of contributions on: <br> <https://github.com/mhinsch/rpathsonpaths/graphs/contributors>

Contributions are welcome via **pull requests**.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

**Maintainer:** Martin Hinsch (<hinsch.martin@gmail.com>)
