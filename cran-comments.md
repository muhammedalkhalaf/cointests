## cointests 1.1.1

CRAN version: 1.1.0.

### Changes

* Documentation only; computations and results are identical to 1.1.0.
* Help page of xtcadfcoint(): the argument simulate is described as the number of Monte Carlo replications used to simulate critical values from independent Gaussian random walks at the estimated break dates.

### Test environments

* Local: R 4.3.3 on Linux (x86_64), R CMD check --as-cran.
* win-builder: R-devel (2026-10-30 r90805 ucrt), Status OK.

### R CMD check results

R CMD check --as-cran on R 4.3.3 (Linux, LANG=C.UTF-8, _R_CHECK_SYSTEM_CLOCK_=0): 0 errors, 0 warnings, 1 NOTE. The NOTE is from the network check: the build environment could not reach orcid.org, so the ORCID iD was reported as possibly invalid; it is the iD in the published 1.1.0.

win-builder, R-devel: Status OK for the same sources with Date 2026-11-01.
