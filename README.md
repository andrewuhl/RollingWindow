# RollingWindow
---

## Intro

The main purpose of this package is to calculate rolling window statistics **fast**. It is aimed at any users who need to calculate rolling statistics
on large data sets, and should be particularly useful for the types of
analysis done in the field of quantitative finance, even though the functions implemented are very general.

## Installing in [R]

```R
# install.packages("devtools") # if not installed

library(devtools)
install_github("andy3rdworld/RollingWindow")
```

## Rolling Window Statistics

This package implements functions which efficiently calculate rolling window statistics in linear time using **single-pass** algorithms. Statistics implemented include:

* beta
* compounding
* covariance
* correlation
* kurtosis
* mean
* mean absolute error
* mean of squares
* mean squared error
* median
* min
* max
* product
* root mean squared error
* skewness
* standard deviation
* sum
* sum product
* sum of squares
* sum of squared errors
* variance
* z score.

_Note_: median is an exception as it is not calculated as a single pass algorithm, although it is still calculated quickly in O(n log window) using efficient algorithms.

Calling functions on rolling windows calculated within loops in [R] can be quite slow, even when proper attempts are made to vectorize operations within the loop. Due to performance issues, several other [R] packages implement commonly used
rolling window functions -- e.g. rolling mean or rolling standard deviation -- using the languages C, Fortran or C++. Both in pure [R] and in compiled language [R] packages, the implementation is usually the naive algorithm; each time a moving window moves forward by one increment, the function is recalulated over all values in the current window.

In this package, rolling window calculations are performed in a single pass, rather than using the naive algorithm. As a rolling window is moved forward, only the oldest value has to be eliminated from the window to calculate the new rolling window statistic. Various known single pass algorithms are used, including a version of Welford's algorithm for
calculating the variance. For minimum and maximum, a double-ended queue (dequeue)-based algorithm is employed. For median, for the rolling window case a circular (ring) buffer data structure is used in its algorithm which achieves O(n log window); for the expanding window case, two balanced priority queues are utilized.

This package's underlying functions are written in C++ and are built on top of the [Rcpp](https://github.com/RcppCore/Rcpp) package. In a simple set of
benchmarks comparing the performance of this package's rolling functions vs. a similar package (also using Rcpp) that
implements a subset of the same functions using the naive algorithm, this package's functions typically complete significantly faster, depending on the window size and function called. For large vector and/or window sizes, it is common to see a 100X+ performance improvement.

Further improvement in performance can be gained by parallelizing the rolling functions, and is left for future work.

## Usage

```R
# General form for single input (vector, matrix, list, etc.)
# result <- RollingFunction(x, window = 10)

# window size
w <- 10

# input vector
x <- rnorm(50)

# only input vector and window size required
avg <- RollingMean(x, window = w)

# std. deviation for population (n)
stddev_sample <- RollingStd(x, window = w, pop = FALSE)

# std. deviation for sample (n - 1)
stddev_pop <- RollingStd(x, window = w, pop = TRUE)

#-------------------------------------------------------------

# General form for functions with two input vectors
# result <- RollingFunction(x, y, window, ...)

# window size
w <- 10

# input vectors
x <- rnorm(50)
y <- rnorm(length(x))

corr    <- RollingCorr(x, y, window = w)
covar   <- RollingCov(x, y, window = w)
beta    <- RollingBeta(x, y, window = w)
sumprod <- RollingSumprod(x, y, window = w)

```

## Benchmarks - Setup

To get a sense of the performance improvement by calculating rolling window
statistics in a single pass, we can run the following in [R].

```R

library(RollingWindow)
library(RcppRoll)
library(microbenchmark)
library(rbenchmark)

# setup
set.seed(10)
n    <- 10000
x    <- rnorm(n); y = rnorm(n)
win  <- 1000
reps <- 10

# relative performance
cols <- c("test", "replications", "elapsed", "relative", "user.self")

```

### Benchmarks - Relative Performance vs. other C++ Implementation
For the first benchmark, let's compare the relative performance of this package's rolling window functions (written in C++) vs. another package's non-single pass implementation (also written in C++).

```R
benchmark(RollingMax(x, win),    roll_max(x, win), replications = reps, columns = cols)

benchmark(RollingMean(x, win),   roll_mean(x, win), replications = reps, columns = cols)

benchmark(RollingMedian(x, win), roll_median(x, win), replications = reps, columns = cols)

benchmark(RollingMin(x, win),    roll_min(x, win), replications = reps, columns = cols)

benchmark(RollingProd(x, win),   roll_prod(x, win), replications = reps, columns = cols)

benchmark(RollingVar(x, win),    roll_var(x, win), replications  = reps, columns = cols)

benchmark(RollingStd(x, win),    roll_sd(x, win), replications    = reps, columns = cols)

```

In the following table, the "relative" column shows relative performance.
The row with **1.0** is the faster of the two functions. For example, in the
standard deviation benchmark, RollingWindow's RollingStd() function
runs 121X faster than RcppRoll's roll_sd() function.

### Max

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingMean | 10 | 0.005 | **1.0** | 0.005 |
| roll_mean | 10 | 0.546 | 109.2 | 0.547 |

### Mean

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingMean | 10 | 0.003 | **1.0** | 0.004 |
| roll_mean | 10 | 0.102 | 34 | 0.101 |

### Median

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingMedian | 10 | 0.018 | **1.0** | 0.019 |
| roll_median | 10 | 6.362 | 353.4 | 6.363 |

### Minimum

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingMin | 10 | 0.005 | **1.0** | 0.005 |
| roll_min | 10 | 0.688 | 137.6 | 0.686 |

### Product

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingProd | 10 | 0.004 | **1.0** | 0.003 |
| roll_prod | 10 | 0.170 | 42.5 | 0.165 |

### Standard Deviation

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingStd | 10 | 0.005 | **1.0** | 0.005 |
| roll_sd | 10 | 0.606 | 121.2 | 0.605 |

### Variance

| test | replications | elapsed | relative | user.self |
| ---- | ------------: | -------: | --------: | ---------: |
| RollingVar | 10 | 0.004 | **1.0** | 0.004 |
| roll_var | 10 | 0.535 | 133.8 | 0.534 |

In all cases, RollingWindow's functions are much faster, between 34X and 353X.

## Benchmarks - Absolute Performance vs. other C++ Implementation

The columns below show the minimum, lower quartile, mean, median, upper quartile and max runtimes, over 10 trials.

```R
Unit: microseconds
          expr        min         lq        mean      median         uq        max neval
    RollingMin    463.788    473.049    714.1995    553.1895    580.332   2350.546    10
      roll_min  68121.511  68135.456  68453.4996  68219.4405  68661.142  69307.602    10
    RollingMax    468.912    485.468    517.9764    525.0100    547.139    572.501    10
      roll_max  54536.383  54557.621  54764.4425  54572.1355  54655.649  56366.473    10
    RollingVar    411.392    428.191    645.7962    490.8550    511.573   2210.802    10
      roll_var  54394.805  54520.076  62992.3953  56243.7960  56426.747 130175.704    10
    RollingStd    515.981    520.592    534.3961    523.9755    530.948    581.469    10
       roll_sd  54169.392  54438.932  55226.7215  54532.8320  56237.989  56771.980    10
    RollingSum    272.602    279.169    480.0387    295.4400    322.925   2136.182    10
      roll_sum  10102.570  10110.305  10163.8197  10150.5820  10188.679  10259.007    10
   RollingMean    325.108    329.620    367.5140    337.5415    425.707    477.601    10
     roll_mean  10102.594  10129.014  10153.1129  10149.5655  10185.238  10200.878    10
   RollingProd    343.545    347.867    397.7944    408.7895    443.519    448.813    10
     roll_prod  16671.556  16675.133  16697.4066  16692.6145  16712.746  16755.501    10
 RollingMedian   1783.405   1795.463   1845.6817   1851.4000   1885.285   1935.361    10
   roll_median 634023.656 635203.964 637219.3584 636503.7700 640487.120 641341.506    10
```
Note that the time unit of the result is _microseconds_. If we look at the median, RollingMedian's mean runtime is 1.8 milliseconds, vs roll_median's mean runtime of 637 milliseconds. While both are fast in the absolute, it is easy to imagine having to run the analysis on several thousand time series with, e.g. daily data. Suppose in a financial application, 2,000 stocks with decades of daily data had to be analyzed. RollingMedian would complete in less than 4 seconds, while roll_median would take 21 minutes. Trying to calculate the same using pure [R] code would be _much_ slower still.

## Benchmarks - All RollingWindow functions

Continuing with the same benchmark setup, the runtimes for all of RollingWindow packages function are see below.

```R
Unit: microseconds
     expr      min       lq      mean    median       uq       max neval
     Beta  564.435  565.453  740.8940  567.8535  570.908  2297.228    10
 Compound  871.097  879.136 3212.4075  914.9325 1098.205 23471.614    10
     Corr  780.402  788.250  803.2185  794.3385  808.005   883.952    10
      Cov  193.533  194.854  203.5932  196.9480  205.787   249.158    10
     Kurt 1420.671 1432.081 1607.6574 1486.9960 1740.712  2090.195    10
      MAE   99.221  100.370 1522.4156  104.1070  109.931 14285.578    10
      Max  467.282  469.833  679.8271  476.0380  556.975  2179.619    10
     Mean  324.000  328.109  374.3722  360.7535  412.956   469.044    10
   Median 1786.301 1790.971 2067.2714 1823.8650 2008.468  3743.916    10
      Min  470.625  472.322  504.8115  483.9675  527.846   591.188    10
       MS  325.019  327.451  340.7884  335.2955  339.075   402.493    10
      MSE   96.127   99.681  117.8423  110.9315  136.579   155.850    10
     Prod  342.187  346.440  364.6057  351.4710  386.365   398.657    10
     RMSE  204.730  208.545  216.3507  217.4755  223.735   229.810    10
     Skew 2105.833 2125.312 2520.8847 2153.6520 2434.727  4798.113    10
       SS  275.331  275.563  289.2251  281.1005  286.441   340.388    10
      SSE   72.808   73.959   79.9912   77.1550   85.095    93.061    10
      Std  518.154  522.499  537.0776  530.7770  539.691   574.066    10
      Sum  273.344  277.292  322.1087  282.8375  351.702   432.871    10
  Sumprod   46.103   48.310   54.2985   49.1690   56.976    80.679    10
      Var  410.086  415.107  457.9157  424.5965  487.923   600.642    10
   Zscore  590.037  596.305  632.3327  605.2265  623.487   841.082    10
```
