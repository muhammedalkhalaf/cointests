## cointests 1.1.1

This is a documentation-only update; the version on CRAN is 1.1.0. There is no numerical change: all statistics, critical values and decisions are identical to 1.1.0.

* Documentation only; no computation changed and all results are identical
  to 1.1.0. In `xtcadfcoint()` the description of the argument `simulate`
  said "bootstrap replications"; it now says that `simulate` is the number
  of Monte Carlo replications used to simulate critical values from
  independent Gaussian random walks at the estimated break dates, as the
  Details section of the help page already stated, and that it is not a
  bootstrap of the data.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
