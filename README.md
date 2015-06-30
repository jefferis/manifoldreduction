# manifoldreduction
## Quick Start

For the impatient ...

```r
# install
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferis/manifoldreduction")

# use
library(manifoldreduction)

# run examples
example("xxx")

# get overview help for package
?manifoldreduction
# help for functions
?xxx

# run tests
library(testthat)
test_package("manifoldreduction")
```

## Installation
Currently there isn't a released version on [CRAN](http://cran.r-project.org/).

### Development version
You can use the **devtools** package to install the development version:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferis/manifoldreduction")
```

Note: Windows users need [Rtools](http://www.murdoch-sutherland.com/Rtools/) and [devtools](http://CRAN.R-project.org/package=devtools) to install this way.
