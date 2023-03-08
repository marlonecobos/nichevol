nichevol: Tools for Ecological Niche Evolution Assessment Considering
Uncertainty
================
Marlon E. Cobos, Hannah L. Owens, and A. Townsend Peterson

- <a href="#package-description" id="toc-package-description">Package
  description</a>
- <a href="#installing-the-package"
  id="toc-installing-the-package">Installing the package</a>
  - <a href="#stable-version" id="toc-stable-version">Stable version</a>
  - <a href="#latest-version" id="toc-latest-version">Latest version</a>
- <a href="#exploring-the-nichevol-package"
  id="toc-exploring-the-nichevol-package">Exploring the nichevol
  package</a>
  - <a href="#setting-a-directory" id="toc-setting-a-directory">Setting a
    directory</a>
  - <a href="#loading-the-package" id="toc-loading-the-package">Loading the
    package</a>
  - <a href="#functions-in-nichevol"
    id="toc-functions-in-nichevol">Functions in nichevol</a>
- <a href="#using-nichevol-with-simple-examples"
  id="toc-using-nichevol-with-simple-examples">Using nichevol with simple
  examples</a>
  - <a href="#packages-needed-for-data-management"
    id="toc-packages-needed-for-data-management">Packages needed for data
    management</a>
  - <a href="#initial-data-example-data"
    id="toc-initial-data-example-data">Initial data (example data)</a>
  - <a href="#preparing-data-for-analyses"
    id="toc-preparing-data-for-analyses">Preparing data for analyses</a>
  - <a href="#ancestral-reconstructions-and-smoothing-of-results"
    id="toc-ancestral-reconstructions-and-smoothing-of-results">Ancestral
    reconstructions and smoothing of results</a>
  - <a href="#representations-of-results"
    id="toc-representations-of-results">Representations of results</a>
- <a href="#references" id="toc-references">References</a>

<br>

<!-- badges: start -->

[![R build
status](https://github.com/marlonecobos/nichevol/workflows/R-CMD-check/badge.svg)](https://github.com/marlonecobos/nichevol/actions)
[![Travis build
status](https://travis-ci.org/marlonecobos/nichevol.svg?branch=master)](https://travis-ci.org/marlonecobos/nichevol)
<!-- badges: end -->

<hr>

## Package description

The **nichevol** R package helps users to perform critical steps in the
process of assessment of species’ ecological niche evolution, with
uncertainty incorporated explicitly in reconstructions. The method
proposed here for ancestral reconstruction of ecological niches
characterizes niches using a bin-based approach that incorporates
uncertainty in estimations. Compared to other existing methods, the
approaches presented here reduce risks of overestimation of amounts or
rates of ecological niche evolution. The main analyses include: initial
exploration of environmental data in occurrence records and accessible
areas, preparation of data for phylogenetic analyses, comparative
phylogenetic analyses, and plotting for interpretation.

<br>

<hr>

## Installing the package

### Stable version

The stable version of **nichevol** is in **CRAN**, and can be installed
using the code below (we are working on this):

``` r
install.packages("nichevol")
```

### Latest version

The most recent version of **nichevol** is available from this GitHub
repository and can be installed using the code below. Please, have in
mind that updates will be done on this version continuously.

Note: Try the code below first… If you have any problem during the
installation, restart your R session, close other sessions you may have
open, and try again. If during the installation you are asked to update
packages, please do it. If any of the packages gives an error, please
install it alone using *install.packages()*, then try re-installing
**nichevol** again. Also, it may be a good idea to update R and RStudio
(if you are using this last).

``` r
# Installing and loading packages
if(!require(remotes)){
    install.packages("remotes")
}
if(!require(nichevol)){
    remotes::install_github("marlonecobos/nichevol")
}
```

<br>

<hr>

## Exploring the nichevol package

### Setting a directory

Some of the main functions of **nichevol** use data that need to be
loaded from a local directory and others produce results that need to be
written to a local directory. Loading data from a local directory and
writing the results outside the R environment helps to avoid problems
related to RAM limitations. That is why setting a working directory is
recommended before starting, as follows:

``` r
directory <- "DRIVE:/YOUR/DIRECTORY" # change the characters accordingly
setwd(directory) 
```

<br>

### Loading the package

Once **nichevol** is installed, you can load the package with the
following line.

``` r
library(nichevol)
```

<br>

### Functions in nichevol

Three main types of functions are included in **nichevol**: (1) ones
that help in preparing data (exploration plots and tables of characters)
for ancestral reconstruction; (2) ones that perform the ancestral
reconstructions (maximum parsimony and maximum likelihood); and (3) some
complementary functions that help in performing post-reconstruction
steps (reconstruction smoothing, and niche and niche evolution
representations). Of course, other helper functions are used in the
package, but they won’t be used as commonly.

A complete list of the functions in the **nichevol** package can be
found in the package documentation. Use the following code to see the
list.

``` r
help(nichevol)
```

<br>

#### Functions for data preparation

These functions are used to explore numerically and graphically the
environments of the areas accessible to the species (**M**) and their
occurrences. They also help users to prepare tables of characters that
represent species’ ecological niches considering used and non-used
conditions, as well as conditions where the use is uncertain. Most of
the functions in this module can consider all species of interest and
multiple environmental variables at the time. For that reason, they read
data from a local directory and have the option to write results to such
directories as well. The functions that work with data from the R
environment are the ones specifically designed to work with multiple
species but only one variable. These last functions do not write results
to local directories. We have intentionally designed some of our
functions to work interacting with local directories to avoid
RAM-related limitations (especially when working with multiple
environmental raster layers at high resolution).

#### Functions for ancestral reconstruction

This module contains functions that help in performing ancestral
reconstruction of species’ ecological niches. These functions use as
inputs the results of the ones from the previous module (tables of
characters) and phylogenetic trees, as in objects of class “phylo” (see
the package **ape**). There are two types of reconstructions available
to date (maximum parsimony and maximum likelihood), but at least one
other type will be included. All these functions use inputs and produce
outputs in the R environment.

#### Functions for post-reconstruction processes

Functions in this module are intended to help with two main processes.
First, one of these functions helps in smoothing results from ancestral
reconstructions. This is necessary to prevent misinterpretations of
results from comparing reconstructed niches of ancestors with niches of
descendants. Other functions help in representing results of previous
analyses. For instance, they help in producing bar-like labels that
represent the niches of species, or they can be used to represent how
niches have evolved across the phylogeny.

<br>

<hr>

## Using nichevol with simple examples

### Packages needed for data management

The following packages are needed for specific tasks. They are used
internally by **nichevol**, and parts of the code for the examples below
will require them. Notice that **nichevol** is already loaded, but these
other packages need to be loaded separately.

``` r
library(terra) # for reading environmental layers and spatial objects
library(ape) # for plotting phylogenetic trees and node labels
library(geiger) # for merging a phylogenetic tree with a table of niche characters
```

### Initial data (example data)

The following lines of code will help to get example data prepared for
demonstrating the usage of **nichevol**. These data are similar to the
ones used in an article in which the methods implemented in **nichevol**
were proposed, illustrated, and explained (see Owens et al. 2020). These
data are included in the package, so their use is straightforward.

``` r
## list of species records
data("occ_list", package = "nichevol")

## list of species accessible areas
m_files <- list.files(system.file("extdata", package = "nichevol"),
                      pattern = "m\\d.gpkg", full.names = TRUE)

m_list <- lapply(m_files, terra::vect)

## raster variable
temp <- rast(system.file("extdata", "temp.tif", package = "nichevol"))

# a simple phylogenetic tree for demonstrations
data("tree", package = "nichevol")
```

<br>

### Preparing data for analyses

Before starting to play with the functions, consider that **nichevol**
allows distinct ways to prepare data, depending on the user’s needs. The
example data downloaded before can be used with the functions designed
to work with multiple variables and all taxa at a time
(`histograms_env`, `stats_evalues`, `bin_tables`, `bin_tables0`).
However, examples with the functions that work with data in the R
environment and only for one variable at a time are shown in detail
here.

#### Exploring data numerically

First check the function documentation:

``` r
help(stats_eval)
```

Now, to run the code,

``` r
stat <- stats_eval(stats = c("mean", "sd", "median", "range", "quantile"), 
                   Ms = m_list, occurrences = occ_list, species = "species",
                   longitude = "x", latitude = "y", variable = temp, 
                   percentage_out = 0)

knitr::kable(stat[[1]], caption = "Table of descriptive statistics of temperature (x 10) in accessible areas for the species of interest.", digits = 2) # %>% kable_styling(font_size = 12)
```

<table>
<caption>
Table of descriptive statistics of temperature (x 10) in accessible
areas for the species of interest.
</caption>
<thead>
<tr>
<th style="text-align:left;">
Species
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
median
</th>
<th style="text-align:right;">
range1
</th>
<th style="text-align:right;">
range2
</th>
<th style="text-align:right;">
quantile.0.
</th>
<th style="text-align:right;">
quantile.25.
</th>
<th style="text-align:right;">
quantile.50.
</th>
<th style="text-align:right;">
quantile.75.
</th>
<th style="text-align:right;">
quantile.100.
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RD 9830
</td>
<td style="text-align:right;">
21.84
</td>
<td style="text-align:right;">
5.54
</td>
<td style="text-align:right;">
23.81
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
27.01
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
21.69
</td>
<td style="text-align:right;">
23.81
</td>
<td style="text-align:right;">
25.35
</td>
<td style="text-align:right;">
27.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 3351
</td>
<td style="text-align:right;">
24.36
</td>
<td style="text-align:right;">
1.65
</td>
<td style="text-align:right;">
24.84
</td>
<td style="text-align:right;">
17.56
</td>
<td style="text-align:right;">
27.01
</td>
<td style="text-align:right;">
17.56
</td>
<td style="text-align:right;">
23.45
</td>
<td style="text-align:right;">
24.84
</td>
<td style="text-align:right;">
25.54
</td>
<td style="text-align:right;">
27.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6933
</td>
<td style="text-align:right;">
15.67
</td>
<td style="text-align:right;">
7.47
</td>
<td style="text-align:right;">
17.53
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
27.16
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
8.49
</td>
<td style="text-align:right;">
17.53
</td>
<td style="text-align:right;">
21.51
</td>
<td style="text-align:right;">
27.16
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 761
</td>
<td style="text-align:right;">
22.28
</td>
<td style="text-align:right;">
5.35
</td>
<td style="text-align:right;">
23.71
</td>
<td style="text-align:right;">
-3.20
</td>
<td style="text-align:right;">
27.01
</td>
<td style="text-align:right;">
-3.20
</td>
<td style="text-align:right;">
22.52
</td>
<td style="text-align:right;">
23.71
</td>
<td style="text-align:right;">
25.21
</td>
<td style="text-align:right;">
27.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6773
</td>
<td style="text-align:right;">
25.00
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
25.14
</td>
<td style="text-align:right;">
16.67
</td>
<td style="text-align:right;">
27.01
</td>
<td style="text-align:right;">
16.67
</td>
<td style="text-align:right;">
24.41
</td>
<td style="text-align:right;">
25.14
</td>
<td style="text-align:right;">
25.74
</td>
<td style="text-align:right;">
27.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 7516
</td>
<td style="text-align:right;">
20.95
</td>
<td style="text-align:right;">
7.68
</td>
<td style="text-align:right;">
25.12
</td>
<td style="text-align:right;">
-3.72
</td>
<td style="text-align:right;">
28.91
</td>
<td style="text-align:right;">
-3.72
</td>
<td style="text-align:right;">
16.40
</td>
<td style="text-align:right;">
25.12
</td>
<td style="text-align:right;">
26.07
</td>
<td style="text-align:right;">
28.91
</td>
</tr>
</tbody>
</table>

``` r

knitr::kable(stat[[2]], caption = "Table of descriptive statistics of temperature (x 10) in occurrences of the species of interest.", digits = 2) #%>% kable_styling(font_size = 12)
```

<table>
<caption>
Table of descriptive statistics of temperature (x 10) in occurrences of
the species of interest.
</caption>
<thead>
<tr>
<th style="text-align:left;">
Species
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
median
</th>
<th style="text-align:right;">
range1
</th>
<th style="text-align:right;">
range2
</th>
<th style="text-align:right;">
quantile.0.
</th>
<th style="text-align:right;">
quantile.25.
</th>
<th style="text-align:right;">
quantile.50.
</th>
<th style="text-align:right;">
quantile.75.
</th>
<th style="text-align:right;">
quantile.100.
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RD 9830
</td>
<td style="text-align:right;">
25.67
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
25.59
</td>
<td style="text-align:right;">
23.79
</td>
<td style="text-align:right;">
27.01
</td>
<td style="text-align:right;">
23.79
</td>
<td style="text-align:right;">
25.36
</td>
<td style="text-align:right;">
25.59
</td>
<td style="text-align:right;">
26.02
</td>
<td style="text-align:right;">
27.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 3351
</td>
<td style="text-align:right;">
25.49
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
25.54
</td>
<td style="text-align:right;">
23.67
</td>
<td style="text-align:right;">
27.00
</td>
<td style="text-align:right;">
23.67
</td>
<td style="text-align:right;">
25.08
</td>
<td style="text-align:right;">
25.54
</td>
<td style="text-align:right;">
25.90
</td>
<td style="text-align:right;">
27.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6933
</td>
<td style="text-align:right;">
25.95
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
25.91
</td>
<td style="text-align:right;">
24.47
</td>
<td style="text-align:right;">
27.11
</td>
<td style="text-align:right;">
24.47
</td>
<td style="text-align:right;">
25.35
</td>
<td style="text-align:right;">
25.91
</td>
<td style="text-align:right;">
26.72
</td>
<td style="text-align:right;">
27.11
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 761
</td>
<td style="text-align:right;">
25.59
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
25.59
</td>
<td style="text-align:right;">
24.02
</td>
<td style="text-align:right;">
26.89
</td>
<td style="text-align:right;">
24.02
</td>
<td style="text-align:right;">
25.21
</td>
<td style="text-align:right;">
25.59
</td>
<td style="text-align:right;">
26.02
</td>
<td style="text-align:right;">
26.89
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6773
</td>
<td style="text-align:right;">
25.34
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
25.41
</td>
<td style="text-align:right;">
23.61
</td>
<td style="text-align:right;">
26.98
</td>
<td style="text-align:right;">
23.61
</td>
<td style="text-align:right;">
24.82
</td>
<td style="text-align:right;">
25.41
</td>
<td style="text-align:right;">
25.83
</td>
<td style="text-align:right;">
26.98
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 7516
</td>
<td style="text-align:right;">
25.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
25.91
</td>
<td style="text-align:right;">
24.09
</td>
<td style="text-align:right;">
27.72
</td>
<td style="text-align:right;">
24.09
</td>
<td style="text-align:right;">
25.31
</td>
<td style="text-align:right;">
25.91
</td>
<td style="text-align:right;">
26.14
</td>
<td style="text-align:right;">
27.72
</td>
</tr>
</tbody>
</table>

To work with multiple variables check the function `stats_evalues`.

<br>

#### Exploring data graphically

First check the help:

``` r
help(hist_evalues)
```

Now, to run the code,

``` r
hist_evalues(M = m_list[[1]], occurrences = occ_list[[1]], species = "species", 
             longitude = "x", latitude = "y", variable = temp,
             CL_lines = c(95, 99), col = c("blue", "red"))
```

![](README_files/figure-gfm/histogram-1.png)<!-- -->

To work with multiple variables check the function `histograms_env`.

<br>

#### Preparing tables of ecological niche characters

First check the help:

``` r
help(bin_table)
```

Now, to run the code,

``` r
bin_tabl <- bin_table(Ms = m_list, occurrences = occ_list, species = "species",
                      longitude = "x", latitude = "y", variable = temp, 
                      percentage_out = 5, n_bins = 20)

knitr::kable(bin_tabl, caption = "Table characters for ecological niches of the species of interest.") #%>% kable_styling(font_size = 12)
```

<table>
<caption>
Table characters for ecological niches of the species of interest.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
3.93 to 5.179
</th>
<th style="text-align:left;">
5.179 to 6.428
</th>
<th style="text-align:left;">
6.428 to 7.677
</th>
<th style="text-align:left;">
7.677 to 8.926
</th>
<th style="text-align:left;">
8.926 to 10.175
</th>
<th style="text-align:left;">
10.175 to 11.424
</th>
<th style="text-align:left;">
11.424 to 12.673
</th>
<th style="text-align:left;">
12.673 to 13.922
</th>
<th style="text-align:left;">
13.922 to 15.171
</th>
<th style="text-align:left;">
15.171 to 16.42
</th>
<th style="text-align:left;">
16.42 to 17.669
</th>
<th style="text-align:left;">
17.669 to 18.918
</th>
<th style="text-align:left;">
18.918 to 20.167
</th>
<th style="text-align:left;">
20.167 to 21.416
</th>
<th style="text-align:left;">
21.416 to 22.665
</th>
<th style="text-align:left;">
22.665 to 23.914
</th>
<th style="text-align:left;">
23.914 to 25.163
</th>
<th style="text-align:left;">
25.163 to 26.412
</th>
<th style="text-align:left;">
26.412 to 27.661
</th>
<th style="text-align:left;">
27.661 to 28.91
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RD 9830
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 3351
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6933
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 761
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6773
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 7516
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
</tr>
</tbody>
</table>

To work with multiple variables check the functions `bin_tables0` and
`bin_tables`.

<br>

### Ancestral reconstructions and smoothing of results

These functions work with one variable at the time, but users can
perform several analyses in a loop if needed. The variable to be
explored here is temperature.

#### Phylogenetic tree and data

With the following code, the phylogenetic tree will be plotted, and its
nodes will be added.

``` r
# exploring tree
tree$tip.label <- rownames(bin_tabl)
plot.phylo(tree, label.offset = 0.04)
nodelabels()
```

![](README_files/figure-gfm/tree_data-1.png)<!-- -->

<br>

#### Maximum parsimony reconstruction

First check the help:

``` r
help(bin_par_rec)
help(smooth_rec)
```

Now, to run the code,

``` r
# tree and data together
tree_data <- treedata(tree, bin_tabl)

# reconstruction
par_rec_table <- bin_par_rec(tree_data)

# smoothing
s_par_rec_table <- smooth_rec(par_rec_table)

# results
knitr::kable(s_par_rec_table, caption = "Table characters for ecological niches of the species of interest and maximum parsimony reconstructions for their ancestors.") #%>% kable_styling(font_size = 12)
```

<table>
<caption>
Table characters for ecological niches of the species of interest and
maximum parsimony reconstructions for their ancestors.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
3.93 to 5.179
</th>
<th style="text-align:left;">
5.179 to 6.428
</th>
<th style="text-align:left;">
6.428 to 7.677
</th>
<th style="text-align:left;">
7.677 to 8.926
</th>
<th style="text-align:left;">
8.926 to 10.175
</th>
<th style="text-align:left;">
10.175 to 11.424
</th>
<th style="text-align:left;">
11.424 to 12.673
</th>
<th style="text-align:left;">
12.673 to 13.922
</th>
<th style="text-align:left;">
13.922 to 15.171
</th>
<th style="text-align:left;">
15.171 to 16.42
</th>
<th style="text-align:left;">
16.42 to 17.669
</th>
<th style="text-align:left;">
17.669 to 18.918
</th>
<th style="text-align:left;">
18.918 to 20.167
</th>
<th style="text-align:left;">
20.167 to 21.416
</th>
<th style="text-align:left;">
21.416 to 22.665
</th>
<th style="text-align:left;">
22.665 to 23.914
</th>
<th style="text-align:left;">
23.914 to 25.163
</th>
<th style="text-align:left;">
25.163 to 26.412
</th>
<th style="text-align:left;">
26.412 to 27.661
</th>
<th style="text-align:left;">
27.661 to 28.91
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RD 9830
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 3351
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6933
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 761
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6773
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 7516
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
</tbody>
</table>

<br>

#### Maximum likelihood reconstruction

First check the help:

``` r
help(bin_ml_rec)
```

Now, to run the code,

``` r
# reconstruction
ml_rec_table <- bin_ml_rec(tree_data)

# smoothing
s_ml_rec_table <- smooth_rec(ml_rec_table)

# results
knitr::kable(s_ml_rec_table, caption = "Table characters for ecological niches of the species of interest and maximum likelihood reconstructions for their ancestors.", digits = 2) #%>% kable_styling(font_size = 12)
```

<table>
<caption>
Table characters for ecological niches of the species of interest and
maximum likelihood reconstructions for their ancestors.
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
3.93 to 5.179
</th>
<th style="text-align:left;">
5.179 to 6.428
</th>
<th style="text-align:left;">
6.428 to 7.677
</th>
<th style="text-align:left;">
7.677 to 8.926
</th>
<th style="text-align:left;">
8.926 to 10.175
</th>
<th style="text-align:left;">
10.175 to 11.424
</th>
<th style="text-align:left;">
11.424 to 12.673
</th>
<th style="text-align:left;">
12.673 to 13.922
</th>
<th style="text-align:left;">
13.922 to 15.171
</th>
<th style="text-align:left;">
15.171 to 16.42
</th>
<th style="text-align:left;">
16.42 to 17.669
</th>
<th style="text-align:left;">
17.669 to 18.918
</th>
<th style="text-align:left;">
18.918 to 20.167
</th>
<th style="text-align:left;">
20.167 to 21.416
</th>
<th style="text-align:left;">
21.416 to 22.665
</th>
<th style="text-align:left;">
22.665 to 23.914
</th>
<th style="text-align:left;">
23.914 to 25.163
</th>
<th style="text-align:left;">
25.163 to 26.412
</th>
<th style="text-align:left;">
26.412 to 27.661
</th>
<th style="text-align:left;">
27.661 to 28.91
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RD 9830
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 3351
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6933
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 761
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 6773
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RD 7516
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
?
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
?
</td>
</tr>
<tr>
<td style="text-align:left;">
LogLik
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
-3.46573590313359
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
-5.49306144388884
</td>
</tr>
<tr>
<td style="text-align:left;">
Rates
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
855.964596564529
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
557.653395948505
</td>
</tr>
<tr>
<td style="text-align:left;">
SE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
1658231.13989086
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
62372.4765130361
</td>
</tr>
</tbody>
</table>

<br>

### Representations of results

#### Ecological niches of species on the phylogeny

``` r
plot.phylo(tree, label.offset = 0.04)
niche_labels(tree, s_par_rec_table, label_type = "tip", height = 0.8, width = 1.5)
```

![](README_files/figure-gfm/tree_niches-1.png)<!-- -->

#### Reconstructed ecological niches of ancestors

``` r
plot.phylo(tree, label.offset = 0.04)
niche_labels(tree, s_par_rec_table, label_type = "tip_node", height = 0.8, width = 1.5)
```

![](README_files/figure-gfm/an_niches-1.png)<!-- -->

#### Evolution of ecological niches in the group

``` r
plot.phylo(tree, label.offset = 0.04)
niche_labels(tree, s_par_rec_table, label_type = "tip", height = 0.8, width = 1.5)
nichevol_labels(tree, s_par_rec_table, height = 0.8, width = 1.5)
```

![](README_files/figure-gfm/niche_evol-1.png)<!-- -->

#### A more informative plot

``` r
par(mfrow = c(1, 2))
plot.phylo(tree, label.offset = 0.04)
niche_labels(tree, s_par_rec_table, label_type = "tip_node", height = 0.8, width = 1.5)
niche_legend(position = "topleft", cex = 0.6)

plot.phylo(tree, label.offset = 0.04)
niche_labels(tree, s_par_rec_table, label_type = "tip", height = 0.8, width = 1.5)
nichevol_labels(tree, s_par_rec_table, height = 0.8, width = 1.5)
nichevol_legend(position = "topleft", cex = 0.6)
```

![](README_files/figure-gfm/niche_evolfin-1.png)<!-- -->

#### Mapping niches and evolution

Evolution occurred between node 9 and the species RD 6933. Let’s map and
see how things look like in geography.

``` r
# preparing layers to represent niches and evolution
niche9 <- map_nichevol(whole_rec_table = s_par_rec_table, variable = temp, 
                       return = "niche", from = "9")
nichesp <- map_nichevol(whole_rec_table = s_par_rec_table, variable = temp, 
                        return = "niche", from = "RD 6933")

evol_8vssp <- map_nichevol(whole_rec_table = s_par_rec_table, variable = temp, 
                           return = "evolution", from = "9", to = "RD 6933")
nichevol_8vssp <- map_nichevol(whole_rec_table = s_par_rec_table, variable = temp, 
                               return = "nichevol", from = "9", to = "RD 6933")

par(mfrow = c(2, 2))
plot(niche9, main = "Niche node 9")
plot(nichesp, main = "Niche RD 6933")

plot(evol_8vssp, main = "Evolution 9 vs RD 6933")
plot(nichevol_8vssp, main = "Niche node 9, evolution to RD 6933")
```

![](README_files/figure-gfm/map_nichevol-1.png)<!-- -->

<br>

<hr>

## References

- Barve, N., V. Barve, A. Jimenez-Valverde, A. Lira-Noriega, S. P.
  Maher, A. T. Peterson, J. Soberón, and F. Villalobos. 2011. The
  crucial role of the accessible area in ecological niche modeling and
  species distribution modeling. Ecological Modelling 222:1810-1819.
- Machado-Stredel, F., M. E. Cobos, and A. T. Peterson. 2021. A
  simulation-based method for selecting calibration areas for ecological
  niche models and species distribution models. Frontiers of
  Biogeography, 13(4):e48814e.
- Owens, H. L., L. P. Campbell, L. L. Dornak, E. E. Saupe, N. Barve, J.
  Soberón, K. Ingenloff, A. Lira-Noriega, C. M. Hensz, C. E. Myers,
  and A. T. Peterson. 2013. Constraints on interpretation of ecological
  niche models by limited environmental ranges on calibration areas.
  Ecological Modelling 263:10-18.
- Owens, H. L., V. Ribeiro, E. E. Saupe, M. E. Cobos, P. A.
  Hosner, J. C. Cooper, A. M. Samy, V. Barve, N. Barve, C. J.
  Muñoz-R, A. T. Peterson. 2020. Acknowledging Uncertainty in
  Evolutionary Reconstructions of Ecological Niches. Ecology and
  Evolution 10(14):6967–6977.
