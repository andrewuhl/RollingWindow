\name{RollingWindow-package}
\alias{RollingWindow}
\docType{package}
\title{
Single pass rolling window calculations
}
\description{
This package implements functions which efficiently calculate rolling window statistics in linear time using single-pass algorithms. Current functions include: beta, compounding, covariance, correlation, kurtosis, mean, mean absolute error, mean of squares, mean squared error, median, min, max, product, root mean squared error, skewness, standard deviation, sum, sum product, sum of squares, sum of squared errors, variance, and z score. Note: median is an exception as it is not calculated as a single pass algorithm, although it is still calculated quickly in O(n log window) using efficient algorithms.
    
Calling functions on rolling windows calculated within loops in [R] can be quite slow, even when proper attempts are made to vectorize operations within the loop. Due to performance issues, several other [R] packages implement commonly used 
rolling window functions -- e.g. rolling mean or rolling standard deviation -- using the languages C, Fortran or C++. Both in pure [R] and in compiled language [R] packages, the implementation is usually the naive algorithm; each time a moving window moves forward by one increment, the function is recalulated over all values in the current window. 

In this package, rolling window calculations are performed in a single pass, rather than using the naive algorithm. As a rolling window is moved forward, generally only the oldest value has to be eliminated from the window to calculate the new rolling window statistic. Various known single pass algorithms are used, including a version of Welford's algorithm for 
calculating the variance. For minimum and maximum, a double-ended queue (dequeue)-based algorithm is employed. For median, for the rolling window case a circular (ring) buffer data structure is used in its algorithm which achieves O(n log window); for the expanding window case, two balanced priority queues are utilized.

This package's underlying functions are written in C++ and are built on top of the Rcpp package. In a simple set of 
benchmarks comparing the performance of this package's rolling functions vs. a similar package (also using Rcpp) that 
implements a subset of the same functions using the naive algorithm, this package's functions typically complete significantly faster, depending on the window size and function called. For large vector and/or window sizes, it is common to see a 100X+ performance improvement.

Further improvement in performance can be gained by parallelizing the rolling functions, and is left for future work.

}

\keyword{ package }
\examples{

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

}

\section{NA Handling}{
For all functions which provide rolling/expanding window statistics in this package, an argument \code{na_rm} can be set to \code{TRUE} to handle \code{NA}'s. 
For functions that accept a single data object \code{x},
 \code{NA}s are handled as follows. For each vector \code{v} in \code{x}:

 1. For all elements \code{i} in \code{v}, if \code{i < window}, then the window statistic for {i} is set to \code{NA}. 
 This is true whether or not the \code{na_rm} argument is set to \code{TRUE}. 
 
 2. For the \emph{rolling} window case, a window statistic is calculated if all elements in the current window are not \code{NA}'s, 
 otherwise the window statistic is set to \code{NA}.
 
 3. For the \emph{expanding} window case, \code{NA}s are ignored.
 
Similarly, for functions accepting two vectors \code{x} and \code{y}:
 
 1. For all elements \code{i} in \code{x} (or equivalently \code{y}, since they are required to be equal-length vectors), 
 if \code{i < window}, then the window statistic for \code{i} is set to \code{NA}. 
 This is true whether or not the \code{na_rm} argument is set to \code{TRUE}. 
 
 2. For the \emph{rolling} window case, a window statistic is calculated if all pairwise elements \code{x_i} and \code{y_i} in the 
 current window are not \code{NA}s, otherwise the window statistic is set to \code{NA}.
 
 3. For the \emph{expanding} window case, \code{NA}'s are ignored if either \code{x_i} or \code{y_i} is an \code{NA}.
 
 It is important to understand what will happen when a rolling statistic is calculated using data that contains \code{NA}s, 
 but when the \code{na_rm} argument is set to \code{FALSE}.
 In this case, for functions like \code{RollingMean}, once the first instance of an \code{NA} is reached, 
 all subsequent elements in the output vector will be \code{NA}. This is because
 most window statistics will include an \code{NA} in an arithmetic or multiplicative operation which will result in an \code{NA}
 for the current window's calculation. However, in other functions like 
 \code{RollingMedian} or \code{RollingMin}, it is not the case that all results including/after the first occurence of an \code{NA} will result in an \code{NA}. 
 This is because arithmetic/multiplicative 
 operations are not performed; instead, comparisons are made between pairs of values which may or may not include an \code{NA}, 
 and any \code{NA}'s resulting from these comparions do not necessarily
 propogate down the result vector.
 
 In short, one should take great care to understand the behavior of \code{NA} handling, and when appropriate, set the \code{na_rm} argument to \code{TRUE}. 
 If the aforementioned \code{NA} handling behavior is not suitable for your application, you might also find the function \code{NaSub} to be useful. \code{NaSub}
 can replace \code{NA}s with either the last non-\code{NA} observation, or a specified (constant) replacement value.
}

\section{Functions arguments x and y}{
If the input object \code{x} contains more than one column of data, then the rolling statistic will be calculated for each column (similar to *apply)
methods in base \code{R}. If a function requires both object \code{x} and a vector \code{y}, the rolling windows statistic will be 
calculated for each column in \code{x} against the vector \code{y}. The dimensions of the return value will be the same as input \code{x}.
}
