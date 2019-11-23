## Test environments
* local windonws 10, R 3.6.1
* local macOS X 10.14.6, R 3.6.1
* ubuntu 16.04.6 LTS (on virtual machine), R 3.6.1
* macOS 10.11 El Capitan (on rhub), R-release
* Windows Server 2008 R2 SP1 32/64 bit (on rhub), R-devel
* Ubuntu 16.04.6 LTS (on travis)


## R CMD check results
There was one ERROR:
* on Ubuntu 16.04.6 LTS (on travis)
When installing package dependencies
  ERROR: configuration failed for package ‘magick’ (libmagick++-dev required)

There were no WARNINGs:

There was one NOTE:

* on Windows Server 2008 R2 SP1 32/64 bit (on rhub), R-devel
Possibly mis-spelled words in DESCRIPTION:
  nichevol (12:14)
  
Action: All words were checked in DESCRIPTION


## Downstream dependencies
There are currently no downstream dependencies for this package.
