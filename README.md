# Cross omic genetic fingerprinting #

![alt text](http://www.molepi.nl/images/logo.png)

This `R` package provides functionality to perform sample relationship
verification using different omic data types. Such as, 450k DNA
methylation data, RNA-sequencing or DNA based genotypes sequencing or
array based imputed.

Checkout the vignettes directory for a few examples.

# Installation #

## Install using the **devtools**-package ##

First install [**devtools**](https://github.com/hadley/devtools). Next
use:

```{r devtools, eval=FALSE}
library(devtools)
install_github("molepi/omicsPrint")
```

Sometimes `install_github` fails with CA cert error. Try running
`httr::set_config(httr::config( ssl_verifypeer = 0L))` before running
`install_github`!

## Install from source using `git/R` ##

Using [git](https://git-scm.com/), e.g, use `git clone` and then build
and install the package from source:

```{r git, engine='bash', eval=FALSE}
git clone git@git.lumc.nl:molepi/omicsPrint.git
R CMD build omicsPrint
R CMD INSTALL omicsPrint_x.y.z.tar.gz
```
Change `_x.y.z.` to the proper version you downloaded!

## Vignettes ##

A few vignettes are provides showing applications of using **omicsPrint** to simulated and publicly available data:

- [omicsPrint](http://htmlpreview.github.io/?https://raw.githubusercontent.com/molepi/omicsPrint/master/inst/doc/omicsPrint.html)

- [GSE102177_DNAm_siblings(http://htmlpreview.github.io/?https://github.com/molepi/omicsPrint/blob/examples/examples/GSE102177_DNAm_siblings.html)

These will be added:
      
GSE52980_DNAm_tissues
GSE73412_DNAm_familymen
GSE39672_DNAm_DNA
GSE100227_DNAm_sisters




    
