## Test environments
* local windonws 10, R 3.6.1
* local macOS X 10.14.6, R 3.6.1
* ubuntu 16.04.6 (on virtual machine), R 3.6.1
* macOS 10.11 El Capitan (on rhub), R-release
* Windows Server 2008 R2 SP1 32/64 bit (on rhub), R-devel


## R CMD check results
There were no ERRORs. 

There were no WARNINGs:

There were four NOTEs:

* on Windows Server 2008 R2 SP1 32/64 bit (on rhub), R-devel
Maintainer: 'Marlon E. Cobos <manubio13@gmail.com>'

License components with restrictions and base license permitting such:
  GPL-3 + file LICENSE

Possibly mis-spelled words in DESCRIPTION:
  nichevol (12:14)

The Description field should not start with the package name,
  'This package' or similar.
  

Actions taken:
* file LICENSE excluded as it was redundat
* all words were checked in DESCRIPTION
* Description field does not start or mention the package name


## Downstream dependencies
There are currently no downstream dependencies for this package.
