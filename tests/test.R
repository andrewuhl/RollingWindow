
library(testthat)
library(RollingWindow)

REPRO = TRUE
TRACE = TRUE
ZAP_DIGITS = 20

veclen = function(x){
  if(is.vector(x))
    length(x)  
  else 
    nrow(x)
}

rand = function(n){
  if(REPRO) set.seed(10)
  rnorm(n)
}

scale_x = function(x) (x - min(x)) / (max(x) - min(x))

#--------------------------------------------
# aliases and/or vectorized methods to match 
# lowercase name of Rolling{func name} function
#--------------------------------------------

# Covariance
Cov = cov

# Variance
Var = var

# Summation
Sum = sum

# Mean
Mean = mean

# Median
Median = median

# Product
Prod = prod

# Minimum
Min = min

# Maximum
Max = max

# standard deviation
Std = sd

# correlation
Corr = cor

# z-score
Zscore = function(x) tail( (x - mean(x)) / sd(x), 1 )

# sum product
Sumprod = function(x, y) as.vector(crossprod(x, y))

# beta
Beta = function(lhs, rhs) cov(lhs, rhs) / var(rhs)

# compounding
Compound = function(x) {
  lenx = veclen(x)
  prod(1+x)^(1/lenx) - 1
}

# mean of squares
MS = function(x) mean(x^2)

# sum of squares
SS = function(x) sum(x^2)

# mean squared error
MSE = function(x, y) mean((x-y)^2)

# sum of squared errors
SSE = function(x, y) sum((x-y)^2)

# mean absolute error
MAE = function(x, y) mean(abs(x-y))

# root mean squared error
RMSE = function(x, y) sqrt(mean((x - y)^2))

# skewness
Skew = function(x, pop = FALSE) {
  M2 = sum((x - mean(x))^2)
  M3 = sum((x - mean(x))^3)
  n = veclen(x)
  res = M3/(M2^(1.5))
  res = res * ifelse(pop, n*sqrt(n-1)/(n-2), sqrt(n))
  return(res)
}

# kurtosis
Kurt = function(x, pop = FALSE) {
  M2 = sum((x - mean(x))^2)
  M4 = sum((x - mean(x))^4)
  n = veclen(x)
  # population excess
  res = n*M4/(M2^2) - 3 
  if(!pop){
    # sample excess
    res = ((n+1)*res +6)*(n-1)/((n-2)*(n-3))
  }
  return(res)
}

# sharpe ratio
Sharpe = function(x, Rf, scale = 12){
  r = prod(1 + x - Rf)^(scale / length(x)) - 1
  res = r / (sd(x) * sqrt(scale))
  return(res)
}

#--------------------------------------------

functions  <- c("Mean", 
                   "Median", 
                   "Prod", 
                   "Min", 
                   "Max", 
                   "Sum",
                   "Std",
                   "Var",
                   "Zscore",
                   "Compound",
                   "MS",
                   "SS",
                   "Skew",
                   "Kurt")

run_tests_x <- function(x, window, functions, trace = FALSE, expanding = FALSE) {
  
  for (i in 1:length(functions)) {
    
    f_name <- functions[i]
    
    if(trace) {
      cat(f_name, "\t\tw =", window, "\texpanding =", expanding, "\n")
    }
    
    f_exp <- get(f_name)
    f_act <- get(paste("Rolling", f_name, sep = ""), envir = asNamespace("RollingWindow"))
    
    exp_ans <- matrix(NA, nrow = nrow(x), ncol = ncol(x), byrow = T)
    if(!is.vector(x)){
      rownames(exp_ans) <- rownames(x)
      colnames(exp_ans) <- colnames(x)  
    }
    for(k in 1:ncol(x)){
      for(w in window:nrow(x)){
        rng <- (w - window + 1):w
        if(expanding) {
          rng <- 1:w
        }
        x_i <- x[rng, k]
        exp_ans[w, k] <- f_exp(x_i)
      }
    }
    exp_ans <- zapsmall(exp_ans, digits = ZAP_DIGITS)
    act_ans <- zapsmall(f_act(x, window, expanding = expanding), digits = ZAP_DIGITS)
    
    test_that(paste(f_name, "matches"), {
      expect_equal(as.vector(act_ans), as.vector(exp_ans))
    })
  }
}

run_tests_xy <- function(x, y, window, functions, trace = FALSE, expanding = FALSE) {
  stopifnot(length(x)==length(y))
  for (i in 1:length(functions)) {
    
    f_name <- functions[i]
    
    if(trace) {
      cat(f_name, "\t\tw =", window, "\texpanding =", expanding, "\n")
    }
    
    f_exp <- get(f_name)
    f_act <- get(paste("Rolling", f_name, sep = ""), envir = asNamespace("RollingWindow"))
    
    exp_ans <- x; exp_ans[] <- NA
    for(k in window:length(x)){
      rng <- (k - window + 1):k
      if(expanding) {
        rng <- 1:k
      }
      x_k <- x[rng]
      y_k <- y[rng]
      exp_ans[k] <- f_exp(x_k, y_k)
    }
    
    exp_ans <- zapsmall(exp_ans, digits = ZAP_DIGITS)
    act_ans <- zapsmall(f_act(x, y, window, expanding = expanding), digits = ZAP_DIGITS)
    
    test_that(paste(f_name, "matches"), {
      expect_equal(as.vector(act_ans), as.vector(exp_ans))
    })
  }
}

# test X vector rolling windows
x <- cbind(rand(50), rand(50))
window <- 5L

run_tests_x(x, window, functions, trace = TRACE)
run_tests_x(x, window, functions, trace = TRACE, expanding = TRUE)

# test when length(x) == window
window <- 50L
run_tests_x(x, window, functions, trace = TRACE)
run_tests_x(x, window, functions, trace = TRACE, expanding = TRUE)

# test functions that should be able to handle window >= 1
window <- 1L
# exclude some functions that can't be calculated for window = 1
funcs = setdiff(functions, c("Std", "Var", "Zscore", "Skew", "Kurt"))
run_tests_x(x, window, funcs)
run_tests_x(x, window, funcs, expanding = TRUE)

# test functions that should be able to handle window >= 2
window <- 2L
# exclude some functions that can't be calculated for window = 2
funcs = setdiff(functions, c("Skew", "Kurt"))
run_tests_x(x, window, funcs, trace = TRACE)
run_tests_x(x, window, funcs, trace = TRACE, expanding = TRUE)

# test copying of row and column names
window <- 10L
rownames(x) <- paste("r", 1:nrow(x), sep = "")
colnames(x) <- paste("c", 1:ncol(x), sep = "")
run_tests_x(x, window, functions, trace = TRACE)
run_tests_x(x, window, functions, trace = TRACE, expanding = TRUE)

# test X, Y vector rolling windows
funcs_xy <- c("Beta", "Cov", "Corr", "Sumprod", "MSE", "SSE", "MAE", "RMSE", "Sharpe")
x <- scale_x(as.matrix(rand(50))) / 100
y <- scale_x(as.matrix(rand(length(x)))) / 100

windows <- c(2, 5, 50)
for(window in windows) {
  run_tests_xy(x, y, window, funcs_xy, trace = TRACE)
  run_tests_xy(x, y, window, funcs_xy, trace = TRACE, expanding = TRUE)
}
