## V 0.1.21 first submission

In this version the following changes were made:

* Niche evolution is labelled as uncertain when ancestral or descendat niche is unknown.
* This change was implemented in functions that add labels and legends to a phylogenetic tree and in the mapping function.

## Test environments

*   Ubuntu 26.04 LTS, R 4.5.3 (local)
*   MacOS 26.7.2, R release (GitHub Actions) 
*   Microsoft Windows Server 2025, R release (GitHub Actions) 
*   Ubuntu 26.04 LTS, R release (GitHub Actions)  
*   Ubuntu 26.04 LTS, R devel (GitHub Actions) 
*   Ubuntu 26.04 LTS, R oldrel-1 (GitHub Actions)  

## R CMD check results

There were no ERRORs.

There were no WARNINGs:

There were no NOTEs:



## Second resubmission 0.1.20
This is a resubmission. In this version I have made the following changes: 

* Added a final point to last sentence of description in DESCRIPTION file, so it is recognized as complete.

**Comments**

* NOTES may be produced because some names of researchers or techniques are considered misspelled words; however, the names are correct.

## Test environments
* Ubuntu 22.04 LTS, R 4.2.2 (local)
* MacOS 12.6.3, R release (GitHub Actions)
* Windows 10.0.20348, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R devel (GitHub Actions)
* Ubuntu 22.04 LTS, R oldrel-1 (GitHub Actions)


## R CMD check results
There were no ERRORs:

There were no WARNINGs:

There were no NOTEs:


## Downstream dependencies
There are currently no downstream dependencies for this package. 

<br>
<hr>

## Resubmission 0.1.20
This is a resubmission. In this version I have made the following changes: 

* Changed package CITATION file to replace old-style personList() or as.personList(). 
* Changed package CITATION file to replace old-style citEntry() with bibentry().
* Replaced URLs by dois in documentation for functions: bin_table, bin_tables0, hist_evalues, histograms_env, stats_eval, and stats_evalues.

## Test environments
* Ubuntu 22.04 LTS, R 4.2.2 (local)
* MacOS 12.6.3, R release (GitHub Actions)
* Windows 10.0.20348, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R devel (GitHub Actions)
* Ubuntu 22.04 LTS, R oldrel-1 (GitHub Actions)


## R CMD check results
There were no ERRORs:

There were no WARNINGs:

There were no NOTEs:


## Downstream dependencies
There are currently no downstream dependencies for this package. 

<br>
<hr>

## V 0.1.20 first submission
In this version I have made the following changes: 

* Replaced functions from the packages sp, raster, and rgdal, with functions from terra.
* Deprecated the argument bin_size and replace it with n_bins to make data processing more efficient in functions: bin_table, bin_tables0, and bin_tables.
* Added the argument verbose to make message suppression easier in functions that needed it: bin_table, bin_tables0, bin_tables, histograms_env, stats_eval, stats_evalues, and map_nichevol (new).
* Added functions to help mapping niche evolution detected based on reconstructions, map_nichevol and helpers.
* Added function to help read prepared data written in a directory: read_bin_table and read_bin_tables.
* Added function to help set uncertainty manually in examples with poorly known species: set_uncertainty.

## Test environments
* Ubuntu 22.04 LTS, R 4.2.2 (local)
* MacOS 12.6.3, R release (GitHub Actions)
* Windows 10.0.20348, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R release (GitHub Actions)
* Ubuntu 22.04 LTS, R devel (GitHub Actions)
* Ubuntu 22.04 LTS, R oldrel-1 (GitHub Actions)


## R CMD check results
There were no ERRORs:

There were no WARNINGs:

There were no NOTEs:


## Downstream dependencies
There are currently no downstream dependencies for this package. 

<br>
<hr>


## Resubmission 0.1.19
This is a resubmission. In this version I have made the following changes:

* Updated the Date field in DESCRIPTION.

**Comments**

* NOTES may be produced because some names of researchers cited in the 
Description file are considered misspelled words; however, the names are correct. 


## Test environments
* local windows 10, R 3.6.2
* local macOS 10.14.6 Mojave, R 3.6.2
* ubuntu 16.04.6 LTS (on travis), R 3.6.2
* macOS 10.11 El Capitan (on rhub), R 3.6.2
* windows server 2008 R2 SP1 32/64 bit (on rhub), R-devel

## R CMD check results
There were no ERRORs:

There were no WARNINGs:

There were no NOTEs:


## Downstream dependencies
There are currently no downstream dependencies for this package. 
