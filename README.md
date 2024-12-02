# madrat based MAgPIE Input Data Library

R package **mrmagpie**, version **1.53.0**

[![CRAN status](https://www.r-pkg.org/badges/version/mrmagpie)](https://cran.r-project.org/package=mrmagpie) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4319612.svg)](https://doi.org/10.5281/zenodo.4319612) [![R build status](https://github.com/pik-piam/mrmagpie/workflows/check/badge.svg)](https://github.com/pik-piam/mrmagpie/actions) [![codecov](https://codecov.io/gh/pik-piam/mrmagpie/branch/master/graph/badge.svg)](https://app.codecov.io/gh/pik-piam/mrmagpie) [![r-universe](https://pik-piam.r-universe.dev/badges/mrmagpie)](https://pik-piam.r-universe.dev/builds)

## Purpose and Functionality

Provides functions for MAgPIE country and cellular input data
    generation.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("mrmagpie")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Kristine Karstens <karstens@pik-potsdam.de>.

## Citation

To cite package **mrmagpie** in publications use:

Karstens K, Dietrich J, Chen D, Windisch M, Alves M, Beier F, Köberle A, v. Jeetze P, Mishra A, Humpenoeder F, Sauer P (2024). _mrmagpie: madrat based MAgPIE Input Data Library_. doi:10.5281/zenodo.4319612 <https://doi.org/10.5281/zenodo.4319612>, R package version 1.53.0, <https://github.com/pik-piam/mrmagpie>.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {mrmagpie: madrat based MAgPIE Input Data Library},
  author = {Kristine Karstens and Jan Philipp Dietrich and David Chen and Michael Windisch and Marcos Alves and Felicitas Beier and Alexandre Köberle and Patrick {v. Jeetze} and Abhijeet Mishra and Florian Humpenoeder and Pascal Sauer},
  year = {2024},
  note = {R package version 1.53.0},
  url = {https://github.com/pik-piam/mrmagpie},
  doi = {10.5281/zenodo.4319612},
}
```
