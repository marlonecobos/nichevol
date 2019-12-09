## Resubmission
This is a resubmission. In this version I have made changes to correct notes 
derived from CRAN check:

There were 4 NOTES:

2 NOTES on r-devel-linux-x86_64-debian-gcc, and 2 on r-devel-windows-ix86+x86_64.
The two notes were the same for both platforms.

* checking CRAN incoming feasibility ... NOTE

  Possibly mis-spelled words in DESCRIPTION:
  phylogenies (16:65)
  
  Package has a VignetteBuilder field but no prebuilt vignette index.

* checking DESCRIPTION meta-information ... NOTE

  Maintainer field differs from that derived from Authors@R
  
  Maintainer: ‘"Marlon E. Cobos" <manubio13@gmail.com>’
  
  Authors@R:  ‘Marlon E. Cobos <manubio13@gmail.com>’
  
  
**Solutions**

The word phylogenies is not mis-spelled.

I excluded the VignetteBuilder field from Description as this package does not
include any vignettes for the moment.

I excluded quotations from the maintainer's name so author and maintainer do not
differ; it was "Marlon E. Cobos", now it is Marlon E. Cobos.

## Test environments
* local windows 10, R 3.6.1
* local macOS X 10.14.6, R 3.6.1
* ubuntu 16.04.6 LTS (on virtual machine), R 3.6.1
* macOS 10.11 El Capitan (on rhub), R-release
* windows server 2008 R2 SP1 32/64 bit (on rhub), R-devel
* ubuntu 16.04.6 LTS (on travis), R 3.6.1


## R CMD check results
There were no ERRORs:

There were no WARNINGs:

There was one NOTE:

* on windows server 2008 R2 SP1 32/64 bit (on rhub), R-devel

Possibly mis-spelled words in DESCRIPTION:
  nichevol (12:14)
  
Action: All words were checked in DESCRIPTION. nichevol is not a grammar issue,
it is the package name.


## Downstream dependencies
There are currently no downstream dependencies for this package.
