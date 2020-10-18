## Description
The recent availability of huge, many-dimensional data sets, like those arising 
from genome-wide association studies (GWAS), provides many opportunities for 
strengthening causal inference. One popular approach is to utilize these 
many-dimensional measurements as instrumental variables (instruments) for 
improving the causal effect estimate between other pairs of variables. 
Unfortunately, searching for proper instruments in a many-dimensional set of 
candidates is a daunting task due to the intractable model space and the fact 
that we cannot directly test which of these candidates are valid. We propose a 
general and efficient causal inference algorithm (MASSIVE) consisting of Model
Assessment and Stochastic Search for Instrumental Variable Estimation. The 
MASSIVE algorithm accounts for model uncertainty by performing Bayesian model 
averaging over the most promising many-dimensional instrumental variable models, 
while at the same time employing weaker assumptions regarding the data 
generating process compared to similar methods.

## Content

The data set contains source code implementing the MASSIVE algorithm, which is 
described in the article titled "[MASSIVE: Tractable and Robust Bayesian Learning 
of Many-Dimensional Instrumental Variable Models](http://proceedings.mlr.press/v124/gabriel-bucur20a.html)" 
by Ioan Gabriel Bucur, Tom Claassen and Tom Heskes. The data set also contains 
simulated data necessary for reproducing the figures in the article as well as 
routines necessary for recreating it. This research is presented in Chapter 5 
of the PhD thesis titled "Being Bayesian about Causal Inference" by
Ioan Gabriel Bucur. The code is written in the R and C++ programming languages.

## Structure

The code is structured on the skeleton of an [R package](https://r-pkgs.org/index.html) 
package as follows:

- The folder `data` contains pre-saved simulated data in .RData format, which we 
use to recreate the figures from the article. The simulated data can also be
reproduced using the `reproduce-data.R` file from the `scripts` folder. 
The simulated data sets are described in `R/data.R`.

- The folder `R` contains the R files that implement the MASSIVE algorithm and the 
routines necessary for reproducing the figures from the article. The main method
is implemented in `R/MASSIVE.R`.

- The folder `man` contains the documentation for the implemented functions.

- The folder `src` contains an efficient Rcpp implementation for computing the
MASSIVE model posterior, gradient, and Hessian matrix.

- The folder `scripts` contains the script `uai-article-figures.R`, which can be 
run from R to produce the figures in the UAI . The script `reproduce-data.R` is 
meant for reproducing the simulated data used in `uai-article-figures.R`.

- The folder `tests` contains a few basic unit tests for the R package.

- The folder `inst` contains additional files used during installation.
  - The subfolder `inst/stan` contains a [Stan](https://mc-stan.org/users/interfaces/stan) 
  file implementing the MASSIVE model, used in `reproduce-data.R` for comparing priors.
  - The subfolder `inst/extdata` contains a few data files that are not in .RData
  format. These contains a list of genetic variants and their associations with
  BMI, the risk of psoriasis, as well as the effect allele frequencies (EAF).
      - `GWAS_BMI.txt` corresponds to Table B (Association of BMI genetic 
    instruments with BMI in UK Biobank and HUNT) in the Supplementary Material of [(Budu-Aggrey et al., 2019)](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002739).
      - `GWAS_psoriasis.txt` corresponds to Table 3 (Association of BMI genetic 
    instruments with psoriasis in UK Biobank and HUNT) in the Supplementary Material of [(Budu-Aggrey et al., 2019)](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002739).
      - `GWAS_EAF.txt` corresponds to Supplementary Table 4 (Significant loci for BMI at P < 1 × 1e-5 for European sex-combined analysis) from the Supplementary Data of [(Locke et al., 2015)](https://www.nature.com/articles/nature14177).
  
- The top folder also contains the following files:
  - `DESCRIPTION` is the file describing the R package
  - `NAMESPACE` is the file specifying the functions provided by the R package
  - `LICENSE.md` is the file containing the GPL-3 license
  
## Prerequisites

In order to install the software, [R](https://cran.r-project.org/) must be 
downloaded and installed. The suggested package [R2BGLiMS](https://github.com/pjnewcombe/R2BGLiMS), 
necessary for running the competing JAM-MR method, also requires the [Java JDK](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) to be installed.
  
## Installation Instructions

Download the software from GitHub with the following command:
`git clone https://github.com/igbucur/MASSIVE.git`. For installing and running the 
MASSIVE R package, several R package are required. These are specified in the 
package `DESCRIPTION` file.

To install the package, open an R instance and run (from the MASSIVE folder):
```
install.packages('devtools') # required package
devtools::install_deps(".", dependencies = TRUE) # install MASSIVE package dependencies
install.packages(".", repos = NULL, type = 'source') # install MASSIVE

library(MASSIVE) # load the package
help(package = "MASSIVE") # see available functions
```

## Licensing

MASSIVE algorithm - Model Assessment and Stochastic Search for Instrumental Variable Estimation

Copyright (C) 2020 Ioan Gabriel Bucur <ioan.gabriel.bucur@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.