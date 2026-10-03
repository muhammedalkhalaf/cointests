# cointests 1.1.1

* Documentation only; no computation changed and all results are identical
  to 1.1.0. In `xtcadfcoint()` the description of the argument `simulate`
  said "bootstrap replications"; it now says that `simulate` is the number
  of Monte Carlo replications used to simulate critical values from
  independent Gaussian random walks at the estimated break dates, as the
  Details section of the help page already stated, and that it is not a
  bootstrap of the data.

# cointests 1.1.0

* `fcoint()` rewritten.
  - FADL is now the ADL regression of Banerjee, Arcabic and Lee (2017) (it was a residual-based test), with lag orders chosen by AIC or BIC on a common sample and the frequency chosen by minimum SSR. Its critical values are those of the Stata module fcoint, which attributes them to Tables 1a and 1b of the paper.
  - FEG follows Yilanci (2019): the ADF regression on the Fourier residuals has no deterministic terms, and the critical values are Table 1 of that paper by number of regressors, frequency and sample size (fixed values were used before).
  - FEG2 adds a constant and the differenced regressors to the ADF regression and reports rho^2; its critical values are not reported because the source tables could not be verified.
  - The Tsong et al. (2016) test follows the TSPDLIB implementation (OLS and DOLS statistics, table of critical values by regressors and frequency). The previous version reported critical values in the wrong order and decided on the wrong tail.
  - FADL, FEG and FEG2 statistics reproduce the Stata module fcoint.
* `xtcadfcoint()` rewritten as a port of the Stata module xtcadfcoint (a translation of the authors' GAUSS code): pooled CCE estimator, break dates by minimum SSR of the defactored (Tb_hat) and detrended (Tb_tilde) residuals, cross-section augmented ADF regressions with impulse dummies, Ng-Perron MAIC and MBIC, and Kim-Perron trimmed statistic for model 5. The previous version averaged unit-level OLS slopes and did not compute the alternative break dates. Results reproduce the Stata module; simulated critical values now use the estimated break dates.
* Corrected the volume and pages of Banerjee and Carrion-i-Silvestre (2025).
* Added tests against Stata reference values.

# cointests 1.0.1

* Corrected the DOI of Banerjee, Arcabic and Lee (2017) to 10.1016/j.econmod.2016.11.004 in DESCRIPTION, R and Rd files.

