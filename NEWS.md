# nichevol 0.1.20 
Fourth release with all tools working and with the following updates:

* Improved algorithms to prepare data for analysis.
* Replaced functions from the packages sp, raster, and rgdal, with functions from terra.
* Deprecated the argument bin_size and replace it with n_bins to make data processing more efficient when creating tables of binary niche characters.
* Improved error and warning messages.
* Added the argument verbose to make message suppression easier in functions that needed it.
* Improved documentation of most functions.
* Added functions to help mapping niche evolution detected based on reconstructions.
* Added function to help read prepared data written in a directory.
* Added function to help set uncertainty manually in examples with poorly known species.


# nichevol 0.1.19 
Third release with all tools working 


# nichevol 0.1.17
Second release after changes suggested by CRAN-member

* Messages are written using `messages()`
* Examples are running while checking
* Examples write info in temporal directory
* Scientific references added to Description

Also, reconstruction smoothing has been improved. 

# nichevol 0.1.16
Initial release after checking
