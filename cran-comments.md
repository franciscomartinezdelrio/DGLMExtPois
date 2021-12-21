## Resubmission
This is a resubmission. In this version I have:

* As suggested by the reviewr, I have removed the redundant "Statistical 
  functions for" from the field Description of the DESCRIPTION file.
  
* Following the reviewer's advice we have included in the Description field
  of the DESCRIPTION file two references to papers describing the methods
  in our package.

* I have removed the \dontrun in the examples in which it was used. I
  used the \dontrun because some examples could take more than 5 secons, but
  it seems that all the examples are fast enough.
  
* As suggested by the reviewr, I have used 
    oldpar <- par(asdf)
    ...
    par(oldpar)
    
  in the examples in which I used the par function.

## Test environments
* local Windows install, R 3.6.1
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There were no ERRORs or WARNINGs. There should be 1 NOTE because this is the first time the package is submitted to CRAN. In two of the platforms in which the package
was tested we found 1 NOTE about a possible mis-spelling in the author's names (Saez and Conde) in DESCRIPTION, but they are correctly spelled.

## Downstream dependencies
The changes made to this package have no effect in downstream dependencies.
