
library(testthat)
library(RollingWindow)

test_that("we act like manually calculated rolling windows", {
  
  TRACE = FALSE
  ZAP_DIGITS = 20
  
  veclen = function(x){
    if(is.vector(x))
      length(x)  
    else 
      nrow(x)
  }
  
  #--------------------------------------------
  # aliases and/or vectorized methods to match 
  # lowercase name of Rolling{func name} function
  #--------------------------------------------
  
  # standard deviation
  std = sd
  
  # correlation
  corr = cor
  
  # z-score
  zscore = function(x) tail( (x - mean(x)) / sd(x), 1 )
  
  # sum product
  sumprod = function(x, y) as.vector(crossprod(x, y))
  
  # beta
  beta = function(x, y) cov(x, y) / var(x)
  
  # compounding
  compound = function(x) {
    lenx = veclen(x)
    prod(1+x)^(1/lenx) - 1
  }
  
  # mean of squares
  ms = function(x) mean(x^2)
  
  # sum of squares
  ss = function(x) sum(x^2)
  
  # mean squared error
  mse = function(x, y) mean((x-y)^2)
  
  # sum of squared errors
  sse = function(x, y) sum((x-y)^2)
  
  # mean absolute error
  mae = function(x, y) mean(abs(x-y))
  
  # root mean squared error
  rmse = function(x, y) sqrt(mean((x - y)^2))
  
  # skewness
  skew = function(x, pop = FALSE) {
    M2 = sum((x - mean(x))^2)
    M3 = sum((x - mean(x))^3)
    n = veclen(x)
    res = M3/(M2^(1.5))
    res = res * ifelse(pop, n*sqrt(n-1)/(n-2), sqrt(n))
    return(res)
  }
  
  # kurtosis
  kurt = function(x, pop = FALSE) {
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
  
  
  #--------------------------------------------
  
  functions_exp <- c("mean", 
                     "median", 
                     "prod", 
                     "min", 
                     "max", 
                     "sum", 
                     "std", 
                     "var",
                     "zscore",
                     "compound",
                     "ms",
                     "ss",
                     "skew",
                     "kurt")
  
  functions_act  <- c("Mean", 
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
  
  run_tests_x <- function(x, window, functions_exp, functions_act, trace = FALSE, expanding = FALSE) {
    stopifnot(length(functions_exp)==length(functions_act))
    for (i in 1:length(functions_exp)) {
      
      f_exp_name <- functions_exp[i]
      f_act_name <- functions_act[i]
      
      if(trace) {
        cat("w =", window, "expanding =", expanding, f_exp_name, "v", f_act_name, "\n")
      }
      
      stopifnot(tolower(f_exp_name) == tolower(f_act_name))
      
      f_exp <- get(f_exp_name)
      f_act <- get(paste("Rolling", f_act_name, sep = ""), envir = asNamespace("RollingWindow"))
      
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
      
      expect_equal(act_ans, exp_ans)
    }
  }
  
  run_tests_xy <- function(x, y, window, functions_exp, functions_act, trace = FALSE, expanding = FALSE) {
    stopifnot(length(functions_exp)==length(functions_act))
    stopifnot(length(x)==length(y))
    for (i in 1:length(functions_exp)) {
      
      f_exp_name <- functions_exp[i]
      f_act_name <- functions_act[i]
      
      if(trace) {
        cat("w =", window, "expanding =", expanding, f_exp_name, "v", f_act_name, "\n")
      }
      
      stopifnot(tolower(f_exp_name) == tolower(f_act_name))
      
      f_exp <- get(f_exp_name)
      f_act <- get(paste("Rolling", f_act_name, sep = ""), envir = asNamespace("RollingWindow"))
      
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
      
      expect_equal(as.vector(act_ans), as.vector(exp_ans))
    }
  }
  
  # test X vector rolling windows
  x <- cbind(rnorm(50), rnorm(50))
  window <- 5L
  
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE)
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE, expanding = TRUE)
  
  # test when length(x) == window
  window <- 50L
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE)
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE, expanding = TRUE)
  
  # test functions that should be able to handle window >= 1
  window <- 1L
  # exclude some functions that can't be calculated for window = 1
  funcsexp = setdiff(functions_exp, c("std", "var", "zscore", "skew", "kurt"))
  funcsact = setdiff(functions_act, c("Std", "Var", "Zscore", "Skew", "Kurt"))
  run_tests_x(x, window, funcsexp, funcsact)
  run_tests_x(x, window, funcsexp, funcsact, expanding = TRUE)
  
  # test functions that should be able to handle window >= 2
  window <- 2L
  # exclude some functions that can't be calculated for window = 2
  funcsexp = setdiff(functions_exp, c("skew", "kurt"))
  funcsact = setdiff(functions_act, c("Skew", "Kurt"))
  run_tests_x(x, window, funcsexp, funcsact, trace = TRACE)
  run_tests_x(x, window, funcsexp, funcsact, trace = TRACE, expanding = TRUE)
  
  # test copying of row and column names
  window <- 10L
  rownames(x) <- paste("r", 1:nrow(x), sep = "")
  colnames(x) <- paste("c", 1:ncol(x), sep = "")
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE)
  run_tests_x(x, window, functions_exp, functions_act, trace = TRACE, expanding = TRUE)
  
  # test X, Y vector rolling windows
  f_xy_exp_name <- c("beta", "cov", "corr", "sumprod", "mse", "sse", "mae", "rmse")
  f_xy_act_name <- c("Beta", "Cov", "Corr", "Sumprod", "MSE", "SSE", "MAE", "RMSE")
  x <- as.matrix(rnorm(50))
  y <- as.matrix(rnorm(length(x)))

  windows <- c(2, 5, 50)
  for(window in windows) {
    run_tests_xy(x, y, window, f_xy_exp_name, f_xy_act_name, trace = TRACE)
    run_tests_xy(x, y, window, f_xy_exp_name, f_xy_act_name, trace = TRACE, expanding = TRUE)
  }
  
})
