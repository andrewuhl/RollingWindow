
library(RollingWindow)
library(RcppRoll)
library(microbenchmark)
library(rbenchmark)

# setup
set.seed(10)
n    <- 1e4
x    <- rnorm(n); y = rnorm(n)
win  <- 1e3
reps <- 10

# relative performance
cols <- c("test", "replications", "elapsed", "relative", "user.self")

benchmark(MAX             = RollingMax(x, win), 
          MAX_RcppRoll    = roll_max(x, win), 
          replications    = reps, columns = cols)

benchmark(MEAN            = RollingMean(x, win), 
          MEAN_RcppRoll   = roll_mean(x, win), 
          replications    = reps, columns = cols)

benchmark(MEDIAN          = RollingMedian(x, win), 
          MEDIAN_RcppRoll = roll_median(x, win), 
          replications    = reps, columns = cols)

benchmark(MIN             = RollingMin(x, win), 
          MIN_RcppRoll    = roll_min(x, win), 
          replications    = reps, columns = cols)

benchmark(PROD            = RollingProd(x, win), 
          PROD_RcppRoll   = roll_prod(x, win), 
          replications    = reps, columns = cols)

benchmark(VAR             = RollingVar(x, win), 
          VAR_RcppRoll    = roll_var(x, win), 
          replications    = reps, columns = cols)

benchmark(STD             = RollingStd(x, win), 
          STD_RcppRoll    = roll_sd(x, win), 
          replications    = reps, columns = cols)

# absolute performance
microbenchmark(MIN           = RollingMin(x, win), 
               MIN_RcppRoll  = roll_min(x, win), 
               
               MAX           = RollingMin(x, win), 
               MAX_RcppRoll  = roll_max(x, win),
               
               VAR           = RollingVar(x, win),
               VAR_RcppRoll  = roll_var(x, win),
               
               STD           = RollingStd(x, win),
               STD_RcppRoll  = roll_sd(x, win),
               
               SUM           = RollingSum(x, win),
               SUM_RcppRoll  = roll_sum(x, win),
               
               MEAN          = RollingMean(x, win),
               MEAN_RcppRoll = roll_mean(x, win), 
               
               PROD          = RollingProd(x, win),
               PROD_RcppRoll = roll_prod(x, win),
               
               MED           = RollingMedian(x, win),
               MED_RcppRoll  = roll_median(x, win), 
               times = reps)

# all of RollingWindow package's functions
microbenchmark(
  Beta     = RollingBeta(x, y, win),
  Compound = RollingCompound(x, win),
  Corr     = RollingCorr(x, y, win),
  Cov      = RollingCov(x, y, win),
  Kurt     = RollingKurt(x, win),
  MAE      = RollingMAE(x, y, win),
  Max      = RollingMax(x, win),
  Mean     = RollingMean(x, win),
  Median   = RollingMedian(x, win),
  Min      = RollingMin(x, win),
  MS       = RollingMS(x, win),
  MSE      = RollingMSE(x, y, win),
  Prod     = RollingProd(x, win),
  RMSE     = RollingRMSE(x, y, win),
  Skew     = RollingSkew(x, win),
  SS       = RollingSS(x, win),
  SSE      = RollingSSE(x, y, win),
  Std      = RollingStd(x, win),
  Sum      = RollingSum(x, win),
  Sumprod  = RollingSumprod(x, y, win),
  Var      = RollingVar(x, win),
  Zscore   = RollingZscore(x, win),
  times = reps)

