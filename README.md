# RollingWindow

## Intro

The purpose of this package is to calculate rolling window and expanding window statistics **fast**. It is aimed at any users who need to calculate rolling statistics
on large data sets, and should be particularly useful for the types of
analysis done in the field of quantitative finance, even though the functions implemented are primarily general-purpose.

## Installing in [R]

```R
# install.packages("devtools") # if not installed

library(devtools)
install_github("andrewuhl/RollingWindow")
```

## Rolling Window Statistics

This package implements functions which efficiently calculate rolling window statistics using mostly **single-pass** algorithms. Statistics implemented include:

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
* sharpe ratio
* skewness
* standard deviation
* sum
* sum product
* sum of squares
* sum of squared errors
* variance
* z score.

_Note_: median and min/max are exceptions as they are not calculated as a single pass algorithm, although they are still calculated quickly in O(n log window) using efficient algorithms.

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
| RollingMax | 10 | 0.005 | **1.0** | 0.005 |
| roll_max | 10 | 0.546 | 109.2 | 0.547 |

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
Note that the time unit of the result is _microseconds_. If we look at the median, RollingMedian's mean runtime is 1.8 milliseconds, vs roll\_median's mean runtime of 637 milliseconds. While both are fast in the absolute, it is easy to imagine having to run the analysis on several thousand time series with, e.g. daily data. Suppose in a financial application, 2,000 stocks with decades of daily data had to be analyzed. RollingMedian would complete in less than 4 seconds, while roll\_median would take 21 minutes. Trying to calculate the same using pure [R] code would be _much_ slower still.

## Benchmarks - All RollingWindow functions

Continuing with the same benchmark setup, the runtimes for all of RollingWindow packages function are seen below.

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

## NA Handling
This package provides support for handling NA's in several ways.

First, the function ```NaSub``` can substitute NA's in a vector by either 
carrying the last non-NA observation forward (similar to the function
```zoo::na.locf```) subject to a maximum gap size between non-NA observations, 
or it can replace NA's with a user-specified replacement value.

### Usage of ```NaSub```
```R
# generate some data with NA's
set.seed(10)
x <- matrix(rnorm(40), nrow = 20)
x[c(2, 5, 7, 9, 14:16), 1] = NA
x[c(1, 2, 15:18), 2] = NA
> print(x)
             [,1]        [,2]
 [1,]  0.01874617          NA
 [2,]          NA          NA
 [3,] -1.37133055 -0.67486594
 [4,] -0.59916772 -2.11906119
 [5,]          NA -1.26519802
 [6,]  0.38979430 -0.37366156
 [7,]          NA -0.68755543
 [8,] -0.36367602 -0.87215883
 [9,]          NA -0.10176101
[10,] -0.25647839 -0.25378053
[11,]  1.10177950 -1.85374045
[12,]  0.75578151 -0.07794607
[13,] -0.23823356  0.96856634
[14,]          NA  0.18492596
[15,]          NA          NA
[16,]          NA          NA
[17,] -0.95494386          NA
[18,] -0.19515038          NA
[19,]  0.92552126 -0.32454401
[20,]  0.48297852 -0.65156299

# replace NA by carrying forward last non-NA observation
y = NaSub(x, last_obs = TRUE, maxgap = 3)
> print(y)
             [,1]        [,2]
 [1,]  0.01874617          NA
 [2,]  0.01874617          NA
 [3,] -1.37133055 -0.67486594
 [4,] -0.59916772 -2.11906119
 [5,] -0.59916772 -1.26519802
 [6,]  0.38979430 -0.37366156
 [7,]  0.38979430 -0.68755543
 [8,] -0.36367602 -0.87215883
 [9,] -0.36367602 -0.10176101
[10,] -0.25647839 -0.25378053
[11,]  1.10177950 -1.85374045
[12,]  0.75578151 -0.07794607
[13,] -0.23823356  0.96856634
[14,] -0.23823356  0.18492596
[15,] -0.23823356          NA
[16,] -0.23823356          NA
[17,] -0.95494386          NA
[18,] -0.19515038          NA
[19,]  0.92552126 -0.32454401
[20,]  0.48297852 -0.65156299

```

In the above example, the first column in ```y``` has had all of
its NA's replaced by the last non-NA observations. However, the 
second column in ```x``` still has NA's remaining in ```y```. For elements in ```x[1:2, 2]```,
this is because the first observation ```x[1, 2]``` was ```NA```, so there
was no non-NA observation to carry forward. For elements in
```x[15:18, 2]```, the NA's were not replaced because the ```maxgap``` argument
was set to 3. This argument limits the maximum gap between non-NA
observations that will be replaced by carrying forward the last non-NA
observation. Since there was a sequence of 4 consecutive NA's, the
NA's were not replaced.

We can call ```NaSub``` with argument ```last_obs = FALSE```. In 
this case, a user-defined value with replace all NA's in lieu
of carrying forward the last non-NA observation, but again subject
to the user-defined ```maxgap``` restriction.

```R
y = NaSub(x, repl = -999, last_obs = FALSE, maxgap = 3)
> print(y)
               [,1]        [,2]
 [1,]    0.01874617          NA
 [2,] -999.00000000          NA
 [3,]   -1.37133055 -0.67486594
 [4,]   -0.59916772 -2.11906119
 [5,] -999.00000000 -1.26519802
 [6,]    0.38979430 -0.37366156
 [7,] -999.00000000 -0.68755543
 [8,]   -0.36367602 -0.87215883
 [9,] -999.00000000 -0.10176101
[10,]   -0.25647839 -0.25378053
[11,]    1.10177950 -1.85374045
[12,]    0.75578151 -0.07794607
[13,]   -0.23823356  0.96856634
[14,] -999.00000000  0.18492596
[15,] -999.00000000          NA
[16,] -999.00000000          NA
[17,]   -0.95494386          NA
[18,]   -0.19515038          NA
[19,]    0.92552126 -0.32454401
[20,]    0.48297852 -0.65156299
```

The ```repl``` argument set the replacement value to ```-999``` ( a number
arbitrarily chosen here for the sake of illustration), and similar
to the first example, column 1 in ```x``` had its NA's replaced by
```-999```, whereas column 2 in ```x``` did not have any NA's replaced due to
the ```maxgap``` argument being set to 3.

If we instead wanted to replace all NA's irrespective of the 
number of consecutive NA's, we need only set the ```maxgap``` to
a larger size. By default, ```maxgap = -1```, which is a special
case indicating that no maximum gap size will be enforced.

```R
y = NaSub(x, repl = -999, last_obs = FALSE, maxgap = -1)
> print(y)
               [,1]          [,2]
 [1,]    0.01874617            NA
 [2,] -999.00000000            NA
 [3,]   -1.37133055   -0.67486594
 [4,]   -0.59916772   -2.11906119
 [5,] -999.00000000   -1.26519802
 [6,]    0.38979430   -0.37366156
 [7,] -999.00000000   -0.68755543
 [8,]   -0.36367602   -0.87215883
 [9,] -999.00000000   -0.10176101
[10,]   -0.25647839   -0.25378053
[11,]    1.10177950   -1.85374045
[12,]    0.75578151   -0.07794607
[13,]   -0.23823356    0.96856634
[14,] -999.00000000    0.18492596
[15,] -999.00000000 -999.00000000
[16,] -999.00000000 -999.00000000
[17,]   -0.95494386 -999.00000000
[18,]   -0.19515038 -999.00000000
[19,]    0.92552126   -0.32454401
[20,]    0.48297852   -0.65156299
```

In this case, all NA's have been replaced except the two leading
NA's in column 2 of ```y```. If we also want to replace leading NA's with
a replacment value, we can set ```leading = TRUE``` and provide
a value for ```repl```. Supposing that we would like leading NA's
to be set to ```-Inf```, we can call ```NaSub``` as follows:

```R
y = NaSub(x, repl = -Inf, last_obs = TRUE, maxgap = -1, leading = TRUE)
             [,1]        [,2]
 [1,]  0.01874617        -Inf
 [2,]  0.01874617        -Inf
 [3,] -1.37133055 -0.67486594
 [4,] -0.59916772 -2.11906119
 [5,] -0.59916772 -1.26519802
 [6,]  0.38979430 -0.37366156
 [7,]  0.38979430 -0.68755543
 [8,] -0.36367602 -0.87215883
 [9,] -0.36367602 -0.10176101
[10,] -0.25647839 -0.25378053
[11,]  1.10177950 -1.85374045
[12,]  0.75578151 -0.07794607
[13,] -0.23823356  0.96856634
[14,] -0.23823356  0.18492596
[15,] -0.23823356  0.18492596
[16,] -0.23823356  0.18492596
[17,] -0.95494386  0.18492596
[18,] -0.19515038  0.18492596
[19,]  0.92552126 -0.32454401
[20,]  0.48297852 -0.65156299
```
In this case, all last non-NA observations were carried forward since
```last_obs = TRUE```, and ```maxgap = -1```. Because ```leading = TRUE```, the leading
values were assigned the ```-Inf``` value set by ```repl = -Inf```. 

### Usage of argument ```na_method``` in ```Rolling*``` functions
```NaSub``` helps replace NA's so that rolling window functions don't have 
to worry about them. So, if ```NaSub``` is able to replace NA's in a helpful way
for your analysis, it makes sense to call ```NaSub``` to transform 
your dataset before calling rolling window functions. 
Nonetheless, the second way to handle NA's is to use the ```na_method``` argument when calling 
```Rolling*``` functions. 

NA handling provided through the ```na_method``` **will not replace NA values**, it only provides
methods of dealing with NA's in calculations (e.g. by ignoring them).

#### na_method = "none"
In this case, no NA handling is performed. If you know that your
data contains no NA's, then ```na_method = "none"``` will be the 
most performant option since it avoids the extra computational cost
and memory allocation to handle NA's.

#### na_method = "ignore"
In this case, NA's are simply ignored. For each column of data in 
dataset ```x```, all NA's are first dropped, then rolling window
statistics are calculated, then a column with the same number of 
observations as the original column containing NA's is written to
with the rolling window statistics, as well as the NA's where they
originally appeared in ```x```. 

#### na_method = "window"
In this case, whenever an NA is first encountered in a vector, a
rolling window statistic will not be calculated again until at least N 
number of consecutive non-NA observations appear, where N is the 
user-specified window length of the rolling window statistic.  

#### Examples of na_method types
The above ```na_method``` types are best illustrated through 
examples. Let's use the rolling window function ```RollingSum```
to illustrate what will be the result of each ```na_method``` type.

```R
x <- matrix(c(1:20, 1:20), ncol = 2)
x <- matrix(c(rep(1, 20), rep(1, 20)), ncol = 2)

> print(x)
      [,1] [,2]
 [1,]    1    1
 [2,]    1    1
 [3,]    1    1
 [4,]    1    1
 [5,]    1    1
 [6,]    1    1
 [7,]    1    1
 [8,]    1    1
 [9,]    1    1
[10,]    1    1
[11,]    1    1
[12,]    1    1
[13,]    1    1
[14,]    1    1
[15,]    1    1
[16,]    1    1
[17,]    1    1
[18,]    1    1
[19,]    1    1
[20,]    1    1

> RollingSum(x, window = 5)
      [,1] [,2]
 [1,]   NA   NA
 [2,]   NA   NA
 [3,]   NA   NA
 [4,]   NA   NA
 [5,]    5    5
 [6,]    5    5
 [7,]    5    5
 [8,]    5    5
 [9,]    5    5
[10,]    5    5
[11,]    5    5
[12,]    5    5
[13,]    5    5
[14,]    5    5
[15,]    5    5
[16,]    5    5
[17,]    5    5
[18,]    5    5
[19,]    5    5
[20,]    5    5
```
As expected, ```RollingSum``` outputs 5 for each element in ```x```
since ```window = 5```. Also, notice that the first 4 (or ```window - 1```)
rows in the result matrix contain NA's. This is by design since we 
only have the required number of observations to calculate the sum 
at row 5 (our specified window length). Suppose now that we add some NA's to 
the two columns in ```x``` and then called ```RollingSum``` once more.

```R
x[11, 1] = NA
x[1, 2] = NA
			
> print(x)
      [,1] [,2]
 [1,]    1   NA
 [2,]    1    1
 [3,]    1    1
 [4,]    1    1
 [5,]    1    1
 [6,]    1    1
 [7,]    1    1
 [8,]    1    1
 [9,]    1    1
[10,]    1    1
[11,]   NA    1
[12,]    1    1
[13,]    1    1
[14,]    1    1
[15,]    1    1
[16,]    1    1
[17,]    1    1
[18,]    1    1
[19,]    1    1
[20,]    1    1

> RollingSum(x, window = 5, method = "none")
      [,1] [,2]
 [1,]   NA   NA
 [2,]   NA   NA
 [3,]   NA   NA
 [4,]   NA   NA
 [5,]    5   NA
 [6,]    5   NA
 [7,]    5   NA
 [8,]    5   NA
 [9,]    5   NA
[10,]    5   NA
[11,]   NA   NA
[12,]   NA   NA
[13,]   NA   NA
[14,]   NA   NA
[15,]   NA   NA
[16,]   NA   NA
[17,]   NA   NA
[18,]   NA   NA
[19,]   NA   NA
[20,]   NA   NA
```

In the first column of output, from the first occurence of NA 
in ```x[, 1]``` at row 11, all subsequent values in the output are NA's. 
This happens because ```na_method = "none"```. Using this method, in the internal C++ calculations each time a 
window is moved forward by one, we add a new value ```x[i, 1]``` to a ```sum``` 
variable at the ```i```th row, then subtract the outgoing value ```x[i - window + 1, 1]``` (i.e. the value that just
left the rolling window ) from the ```sum``` variable to calculate 
the rolling sum for the current observation. Since we've added the new value of NA at 
```x[11, 1]``` to ```sum```, it will from that index forward be equal 
to NA, which is exactly the result we would expect if we implemented a rolling
sum function in R using the same single pass summation algorithm 
described above that doesn't account for NA's.

At this point, we have two ways of avoiding NA's propogating down
the columns of the output. We can set ```na_method = "window"``` or
```na_method = "ignore"```. Let's try both types and see what 
happens.

```R
> RollingSum(x, window = 5, na_method = "window")
      [,1] [,2]
 [1,]   NA   NA
 [2,]   NA   NA
 [3,]   NA   NA
 [4,]   NA   NA
 [5,]    5   NA
 [6,]    5    5
 [7,]    5    5
 [8,]    5    5
 [9,]    5    5
[10,]    5    5
[11,]   NA    5
[12,]   NA    5
[13,]   NA    5
[14,]   NA    5
[15,]   NA    5
[16,]    5    5
[17,]    5    5
[18,]    5    5
[19,]    5    5
[20,]    5    5
```

With ```na_method = "window"```, we get more intuitive results. Once the
first NA occurs in ```x[11, 1]```, the next statistic isn't outputted
until 5 consecutive non-NA values appear (5 being the length of
the window we specified), which first occurs at ```x[16, 1]```.

However, we also have one more option, which is setting 
```na_method = "ignore"```. Let's see what happens in this case.

```R
> RollingSum(x, window = 5, na_method = "ignore")
     [,1] [,2]
 [1,]   NA   NA
 [2,]   NA   NA
 [3,]   NA   NA
 [4,]   NA   NA
 [5,]    5   NA
 [6,]    5    5
 [7,]    5    5
 [8,]    5    5
 [9,]    5    5
[10,]    5    5
[11,]   NA    5
[12,]    5    5
[13,]    5    5
[14,]    5    5
[15,]    5    5
[16,]    5    5
[17,]    5    5
[18,]    5    5
[19,]    5    5
[20,]    5    5
```
Here, NA's only appear where they occur in the original dataset ```x```.
This is because ```na_method = "ignore"``` completely ignores NA's, 
that is, NA's are simply skipped when moving the window forward. The 
sums will still be calculated using exactly 5 data points (the size
of the window). If we look at the output at ```[12, 1]```, the sum
was calculated as ```x[12, 1] + x[10, 1]  + x[9, 1] + x[8, 1] + x[7, 1]```.
Notice that ```x[11, 1]``` was not included, but 5 observations were
included, thus, the window size is always respected. This can be seen
more clearly if we change the values of x to an increasing sequence
of integers.

```R
x <- matrix(c(1:20, 1:20), ncol = 2)
x[11, 1] = NA
x[1, 2] = NA

> print(x)
      [,1] [,2]
 [1,]    1   NA
 [2,]    2    2
 [3,]    3    3
 [4,]    4    4
 [5,]    5    5
 [6,]    6    6
 [7,]    7    7
 [8,]    8    8
 [9,]    9    9
[10,]   10   10
[11,]   NA   11
[12,]   12   12
[13,]   13   13
[14,]   14   14
[15,]   15   15
[16,]   16   16
[17,]   17   17
[18,]   18   18
[19,]   19   19
[20,]   20   20

> RollingSum(x, window = 5, na_method = "ignore")
      [,1] [,2]
 [1,]   NA   NA
 [2,]   NA   NA
 [3,]   NA   NA
 [4,]   NA   NA
 [5,]   15   NA
 [6,]   20   20
 [7,]   25   25
 [8,]   30   30
 [9,]   35   35
[10,]   40   40
[11,]   NA   45
[12,]   46   50
[13,]   52   55
[14,]   58   60
[15,]   64   65
[16,]   70   70
[17,]   75   75
[18,]   80   80
[19,]   85   85
[20,]   90   90
```

Until ```[11, 1]```, the output increased by 5 at each row, but
from ```[12:15, 1]``` -- since the sequence is disrupted (i.e. NA is 
skipped) -- the sequence increases by 6 until there are 5 consecutive
non-NA values, at which point the normal output sequence of 
increasing by 5 occurs at ```[16, 1]```.

#### Expanding window NA handling

Until now, we have only discussed the **rolling window** case
 (```expanding = FALSE```) and not the **expanding window** case 
 (```expanding = TRUE```). For the expanding case, if ```na_method = "none"```,
 there is no support for NA handling similar to the rolling window case.
 For the cases ```na_method = "ignore"``` and ```na_method = "window"```,
 both will use the same logic as ```na_method = "ignore"``` as in 
 the rolling window case.