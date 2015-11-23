/*
 Based on http://opensource.org/licenses/MIT

Copyright (c) 2015, Andrew Uhl

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <sstream>
#include <string>
#include <deque>
#include <queue>
#include <cmath>

#include <Rcpp.h>
#include "rolling_median.h"
#include "expanding_median.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

/*
* TODO:
* 
*  CACHE LOCALITY
*  ->  any way to not reach back in array to improve cache locality?
*  
*  WINDOW SIZE OPTIMIZATION
*  ->  for which window sizes is it optimal to use a naive algorithm vs. single pass?
*  
*  PARALLELIZE
*  ->  add parallel support for:
*        (1) single vector calculations with multi-column matrices, and/or
*        (2) all calculations with vector length > N (N to be determined)
*/

// for one-vector calculations {std. deviation, variance, z-score} // 
enum CalcTypeVar {STDDEV, VAR, ZSCORE};

// for two-vector calculations {covariance, correlation}
enum CalcTypeCov {BETA, COV, CORR};

// for one-vector calculations of higher moments {skewness and kurtosis}
enum CalcTypeMom {SKEW, KURT};

// for squared error calculations {mean of sq err, root mean sq err, sum of sq err}
enum CalcTypeSqErr {MEAN, RMSE, SUM};

// basic error constants
enum ErrorType {XY_LEN_NE, X_LEN_EQ_0, X_LEN_LT_WIN, X_LEN_EQ_1, WIN_LEN_EQ_0, WIN_LEN_EQ_1, WIN_LEN_EQ_N};

// print errors to R console
void print_error(ErrorType etype, int minwin = 1) {
  std::string errmsg = "";
  switch(etype){
  case XY_LEN_NE:
    errmsg = "vector lengths must be equal";
    break;
  case X_LEN_LT_WIN:
    errmsg = "vector length(s) must be <= window size";
    break;
  case X_LEN_EQ_0:
    errmsg = "vector length(s) must be > 0";
    break;
  case X_LEN_EQ_1:
    errmsg = "vector length(s) must be > 1";
    break;
  case WIN_LEN_EQ_0:
    errmsg = "window size must be > 0";
    break;
  case WIN_LEN_EQ_1:
    errmsg = "window size must be > 1";
    break;
  case WIN_LEN_EQ_N:
    std::stringstream sstm;
    sstm << "window size must be > " << minwin;
    errmsg = sstm.str();
    break;
  }
  if (errmsg == "")
    errmsg = "unsepcified error";
  Rcpp::stop(errmsg.c_str());
}

// converts SEXP to NumericMatrix
NumericMatrix as_matrix(const SEXP& x) {
  NumericMatrix y = internal::convert_using_rfunction(x, "as.matrix");
  return y;
}

// assign rownames and colnames from source to target
void setnames(NumericMatrix& source, NumericMatrix& target) {
  target.attr("dimnames") = clone(List(source.attr("dimnames")));
}

// check x for cases when length(x) must be >= n
void check_xn(const NumericVector& x, int window, int minwin){
  if (x.length() < 2)
    print_error(X_LEN_EQ_1);
  
  if (window < minwin)
    print_error(WIN_LEN_EQ_N, minwin);
  
  if (x.length() < window)
    print_error(X_LEN_LT_WIN);
}

// check x for cases when length(x) must be >= 2
void check_x(const NumericVector& x, int window){
  if (x.length() < 2)
    print_error(X_LEN_EQ_1);
  
  if (window < 2)
    print_error(WIN_LEN_EQ_1);
  
  if (x.length() < window)
    print_error(X_LEN_LT_WIN);
}

// check x for cases when length(x) must be >= 1
void check_x1(const NumericVector& x, int window){
  if (x.length() == 0)
    print_error(X_LEN_EQ_0);
  
  if (window == 0)
    print_error(WIN_LEN_EQ_0);
  
  if (x.length() < window)
    print_error(X_LEN_LT_WIN);
}

// check x and y for two-vector calculations
void check_xy(const NumericVector& x, const NumericVector& y, int window){
  if (x.length() != y.length())
    print_error(XY_LEN_NE);
  
  if (window < 2)
    print_error(WIN_LEN_EQ_1);
  
  if (x.length() < 2)
    print_error(X_LEN_EQ_1);
  
  if (x.length() < window)
    print_error(X_LEN_LT_WIN);
}

// check x and y for two-vector calculations
void check_xy1(const NumericVector& x, const NumericVector& y, int window){
  if (x.length() != y.length())
    print_error(XY_LEN_NE);
  
  if (window == 0)
    print_error(WIN_LEN_EQ_0);
  
  if (x.length() < window)
    print_error(X_LEN_LT_WIN);
}

/* ------------------------------
* PRIVATE
* 
* Vector-based calculations
* ------------------------------
*/

// calculates rolling compound (useful for financial asset returns)
NumericVector rolling_compound(const NumericVector& x, int window, long double scale, bool expanding) {
  
  check_x1(x, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double prod = 1;
  for (int i = 0; i < n; ++i) {
    prod *= 1 + x[i];
    if (i >= window - 1){
      if (expanding) {
        rollx[i] = pow(prod, scale / (long double) (i + 1)) - 1;
      } else {
        rollx[i] = pow(prod, scale / (long double) window) - 1;
        prod /= (1 + x[i - w1]);  
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window for one of {covariance, correlation, beta}
NumericVector rolling_covcorrbeta(const NumericVector& x, const NumericVector& y, 
                                  int window, bool pop, bool expanding, CalcTypeCov ctype) {
  
  if (pop) {
    check_xy1(x, y, window); 
  } else {
    check_xy(x, y, window);  
  }
  
  int W = window;
  int n_xy = x.length();
  int w = 0;
  int pop_n = pop ? W : W - 1;
  long double avg_x = 0, sumsq_x = 0, delta_x = 0, var_x = 0, sd_x = 0, sum_x = 0;
  long double avg_y = 0, sumsq_y = 0, delta_y = 0, var_y = 0, sd_y = 0, sum_y = 0;
  long double sum_xy = 0;
  long double cov  = 0;
  
  NumericVector rollxy(n_xy);
  
  if (expanding) {
    for (int i = 0; i < n_xy; ++i) {
      ++w;
      delta_x = x[i] - avg_x;
      avg_x += delta_x / w;
      sumsq_x += delta_x * (x[i] - avg_x);
      sum_x += x[i];
      
      delta_y = y[i] - avg_y;
      avg_y += delta_y / w;
      sumsq_y += delta_y * (y[i] - avg_y);
      sum_y += y[i];
      
      sum_xy += x[i] * y[i];
      
      pop_n = (pop ? w : w - 1);
      var_x = sumsq_x / pop_n;
      var_y = sumsq_y / pop_n;
      if(ctype == CORR) {
        sd_x = sqrt(var_x);
        sd_y = sqrt(var_y);  
      }
      cov = (sum_xy - sum_x * sum_y / w) / pop_n;
      
      if (i >= window - 1) {
        switch (ctype) {
        case BETA : rollxy[i] = cov / var_x; break;
        case CORR : rollxy[i] = cov / (sd_x * sd_y); break;
        case COV  : rollxy[i] = cov; break;
        }  
      } else {
        rollxy[i] = NA_REAL;
      }
    }
    return rollxy;
  }
  
  for (int i = 0; i < W; ++i) {
    ++w;
    delta_x = x[i] - avg_x;
    avg_x += delta_x / w;
    sumsq_x += delta_x * (x[i] - avg_x);
    sum_x += x[i];
    
    delta_y = y[i] - avg_y;
    avg_y += delta_y / w;
    sumsq_y += delta_y * (y[i] - avg_y);
    sum_y += y[i];
    
    sum_xy += x[i] * y[i];
    
    rollxy[i] = NA_REAL;
  }
  
  var_x = sumsq_x / pop_n;
  var_y = sumsq_y / pop_n;
  if(ctype == CORR) {
    sd_x = sqrt(var_x);
    sd_y = sqrt(var_y);  
  }
  cov = (sum_xy - sum_x * sum_y / W) / pop_n;
  
  switch (ctype) {
  case BETA : rollxy[W - 1] = cov / var_x; break;
  case CORR : rollxy[W - 1] = cov / (sd_x * sd_y); break;
  case COV  : rollxy[W - 1] = cov; break;
  }
  
  // std dev terms
  long double xi_old = x[0], xi = 0, avg_old_x = 0;
  long double yi_old = y[0], yi = 0, avg_old_y = 0;
  
  for (int i = W; i < n_xy; ++i) {
    
    // std dev of x
    xi = x[i];
    if (ctype == BETA || ctype == CORR) {
      avg_old_x = avg_x;
      avg_x = avg_old_x + (xi - xi_old) / W;
      var_x += (xi - xi_old)*(xi - avg_x + xi_old - avg_old_x) / pop_n;
      if(ctype == CORR)
        sd_x = sqrt(var_x);
    }
    
    // std dev of y
    yi = y[i];
    if (ctype == BETA || ctype == CORR) {
      avg_old_y = avg_y;
      avg_y = avg_old_y + (yi - yi_old) / W;
      var_y += (yi - yi_old)*(yi - avg_y + yi_old - avg_old_y) / pop_n;
      if(ctype == CORR)
        sd_y = sqrt(var_y);
    }
    
    // cov of x,y
    sum_xy += xi * yi - xi_old * yi_old;
    sum_x += xi - xi_old;
    sum_y += yi - yi_old;
    cov = (sum_xy - sum_x * sum_y / W) / pop_n;
    
    switch (ctype) {
    case BETA : rollxy[i] = cov / var_x; break;
    case CORR : rollxy[i] = cov / (sd_x * sd_y); break;
    case COV  : rollxy[i] = cov; break;
    }
    
    xi_old = x[i - W + 1];
    yi_old = y[i - W + 1];
  }
  return rollxy;
}

// calculates rolling mean
NumericVector rolling_mean(const NumericVector& x, int window, bool expanding) {
  
  check_x1(x, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += x[i];
    if (i >= window - 1){
      if (expanding) {
        rollx[i] = sum / (i + 1);
      } else {
        rollx[i] = sum / window;
        sum -= x[i - w1];  
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling mean absolute error
NumericVector rolling_meanabserr(const NumericVector& x, const NumericVector& y, int window, 
                                 bool expanding) {
  
  check_xy1(x, y, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sumabserr = 0;
  for (int i = 0; i < n; ++i) {
    sumabserr += std::abs(x[i] - y[i]);
    if (i >= window - 1) {
      if (expanding) {
        rollx[i] = sumabserr / (i + 1);
      } else {
        rollx[i] = sumabserr / window;  
        sumabserr -= std::abs(x[i-w1] - y[i - w1]);
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calcuates rolling median
NumericVector rolling_med(const NumericVector& x, int window, bool expanding) {
  
  check_x1(x, window);
  
  int n = x.length();
  NumericVector rollx(n);
  
  if (expanding) {
    
    ExpandingMedian* median = new ExpandingMedian();
    
    for (int i = 0; i < n; ++i) {
      median->Insert(x[i]);
      if (i >= window - 1) {
        rollx[i] = median->GetMedian();  
      } else {
        rollx[i] = NA_REAL;
      }
    }
    delete median;
    
  } else {
    
    Mediator* m = MediatorNew(window);
    
    for (int i = 0; i < n; ++i) {
      MediatorInsert(m, x[i]);
      if (i >= window - 1) {
        rollx[i] = MediatorMedian(m);
      } else {
        rollx[i] = NA_REAL;
      }
    }
    free(m); // malloc inside Mediator  
    
  }
  
  return rollx;
}

// calculates rolling window for one of {minimum, maximum}
NumericVector rolling_minmax(const NumericVector& x, int window, bool min, bool expanding) {
  
  check_x1(x, window);
  
  int n  = x.length();
  NumericVector rollx(n);
  
  if (expanding) {
    long double minmax = x[0];
    for (int i = 0; i < n; ++i) {
      if (min) {
        if (x[i] < minmax) {
          minmax = x[i];
        }
      } else {
        if (x[i] > minmax) {
          minmax = x[i];
        }
      }
      if (i >= window - 1) {
        rollx[i] = minmax;  
      } else {
        rollx[i] = NA_REAL;
      }
    }
    return rollx;
  }
  
  std::deque< std::pair<long double, int> > deck;
  for (int i = 0; i < x.size(); ++i) {
    if(min) {
      while (!deck.empty() && deck.back().first >= x[i])
        deck.pop_back();
    } else {
      while (!deck.empty() && deck.back().first <= x[i])
        deck.pop_back();  
    }
    deck.push_back(std::make_pair(x[i], i));
    
    while(deck.front().second <= i - window)
      deck.pop_front();
    
    long double min = deck.front().first;
    if (i < window - 1) {
      rollx[i] = NA_REAL;
    } else {
      rollx[i] = min;
    }     
  }
  return rollx;
}

// calculates rolling product
NumericVector rolling_prod(const NumericVector& x, int window, bool expanding) {
  
  check_x1(x, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double prod = 1;
  for (int i = 0; i < n; ++i) {
    prod *= x[i];
    if (i >= window - 1) {
      rollx[i] = prod;
      if (!expanding)
        prod /= x[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window for one of {skewness, kurtosis}
NumericVector rolling_skewkurt(const NumericVector& x, int window, bool pop, CalcTypeMom ctype, bool expanding) {
  
  if (ctype == SKEW) {
    if (pop) {
      check_xn(x, window, 3);
    } else {
      check_xn(x, window, 2);
    }
  } else if (ctype == KURT) {
    check_xn(x, window, 4);
  }
  
  // ----------------------------
  // expanding window calculation  
  // ----------------------------
  
  int W = window;
  int n_x = x.length();
  NumericVector rollx(n_x);
  
  if (expanding) {
    
    int n = 0;
    long double delta = 0, delta_n = 0, delta_n2 = 0, term1 = 0;
    long double M1 = 0, M2 = 0, M3 = 0, M4 = 0;
    long double pop_skew = 0;
    
    for (int i = 0; i < n_x; ++i) {
      
      // ----------- BEGIN ATTRIBUTION --------------------------
      // this block borrowed from Jonh D. Cook, MIT license
      //      http://www.johndcook.com/blog/skewness_kurtosis/
      // --------------------------------------------------------
      ++n;
      delta = x[i] - M1;
      delta_n = delta / n;
      delta_n2 = delta_n * delta_n;
      term1 = delta * delta_n * (n - 1);
      M1 += delta_n;
      if (ctype == KURT) {
        M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;  
      }
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
      // ----------- END ATTRIBUTION ----------------------------
      
      if (i < W - 1) {
        rollx[i] = NA_REAL;
      } else {
        if (ctype == SKEW) {
          pop_skew = pop ? double(n)*sqrt(double(n)-1)/(double(n)-2): sqrt(double(n));
          rollx[i] = pop_skew * M3 / pow(M2, 1.5);
        } else if (ctype == KURT) {
          rollx[i] = n * M4 / (M2*M2) - 3;
          if(!pop){
            rollx[i] = ((n+1)*rollx[i] +6)*(n-1)/((n-2)*(n-3));
          }
        }  
      }
    }
    return rollx;
  } else {
    
    long double pop_skew = 0;
    
    // ----------------------------
    // rolling window calculation  
    // ----------------------------
    
    if (ctype == SKEW) {
      pop_skew = pop ? double(W)*sqrt(double(W)-1)/(double(W)-2): sqrt(double(W));
    }
    
    // 2nd moment (M2)
    int n = 0;
    long double avg = 0, delta = 0, M2 = 0;
    long double xi_old = x[0], xi = 0, avg_old = 0;
    
    // 3rd moment (M3)
    long double sumx = 0, sumx2 = 0, sumx3 = 0, x2 = 0, x2_old;
    long double M3 = 0;
    long double c1_M3 = 3.0 / double(W), c2_M3 = 2.0 / pow(double(W), 2);
    int w1 = W - 1;
    
    // 4th moment (M4)
    long double sumx4 = 0;
    long double M4 = 0;
    long double c1_M4 = 6.0 / pow(double(W), 2), c2_M4 = -4.0 / double(W), c3_M4 = -3.0 / pow(double(W), 3);
    
    for (int i = 0; i < n_x; ++i) {
      xi     = x[i];
      x2     = xi*xi;
      sumx  += xi;
      sumx2 += x2;
      sumx3 += x2*xi;
      sumx4 += x2*x2;
      
      if (i < W) {
        ++n;
        delta = xi - avg;
        avg += delta / n;
        M2 += delta * (xi - avg);
        rollx[i] = NA_REAL;
      } else {
        avg_old = avg;
        avg = avg_old + (xi - xi_old) / W;
        M2 += (xi - xi_old)*(xi - avg + xi_old - avg_old);
      }
      if (i >= W - 1) {
        xi_old   = x[i - w1];
        
        if (ctype == SKEW) {
          M3 = sumx3 - c1_M3*sumx*sumx2 + c2_M3*pow(sumx, 3);
          rollx[i] = pop_skew * M3 / pow(M2, 1.5);  
        } else if (ctype == KURT) {        
          M4 = sumx4 + sumx*(c1_M4*sumx*sumx2 + c2_M4*sumx3 + c3_M4*pow(sumx, 3));
          rollx[i] = W*M4 / (M2*M2) - 3;
          if (!pop) {
            rollx[i] = ((W+1)*rollx[i] + 6) * (W-1)/((W-2)*(W-3));
          }
        }
        
        x2_old   = xi_old*xi_old;
        sumx    -= xi_old;
        sumx2   -= x2_old;
        sumx3   -= x2_old*xi_old;
        sumx4   -= x2_old*x2_old;
      }
    }
    return rollx;
    
  }
}

// calculates rolling sum
NumericVector rolling_sum(const NumericVector& x, int window, bool expanding) {
  
  check_x1(x, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += x[i];
    if (i >= window - 1) {
      rollx[i] = sum;
      if (!expanding)
        sum -= x[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling sum product
NumericVector rolling_sumprod(const NumericVector& x, const NumericVector& y, int window, bool expanding) {
  
  check_xy1(x, y, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sumprod = 0;
  for (int i = 0; i < n; ++i) {
    sumprod += x[i] * y[i];
    if (i >= window - 1) {
      rollx[i] = sumprod;
      if (!expanding)
        sumprod -= x[i - w1] * y[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates one of rolling {mean of squares, sum of squares}
NumericVector rolling_sumsq(const NumericVector& x, int window, bool expanding, CalcTypeSqErr ctype) {
  
  check_x1(x, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sumsq = 0;
  for (int i = 0; i < n; ++i) {
    sumsq += x[i] * x[i];
    if (i >= window - 1) {
      rollx[i] = sumsq; 
      if (ctype == MEAN) {
        rollx[i] /= expanding ? i + 1 : window;
      }
      if (!expanding) {
        sumsq -= x[i - w1] * x[i - w1];
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling sum of squared errors
NumericVector rolling_sqerr(const NumericVector& x, const NumericVector& y, int window, 
                            bool expanding, CalcTypeSqErr ctype) {
  
  check_xy1(x, y, window);
  
  int n  = x.length(), w1 = window - 1;
  NumericVector rollx(n);
  long double sumsqerr = 0;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    sumsqerr += x[i]*x[i] - 2*x[i]*y[i] + y[i]*y[i];
    if (i >= window - 1) {
      rollx[i] = sumsqerr;
      if (ctype == MEAN || ctype == RMSE) {
        rollx[i] /= expanding ? i + 1 : window;
        if (ctype == RMSE) {
          rollx[i] = sqrt(rollx[i]);
        }
      }
      if (!expanding) {
        k = i - w1;
        sumsqerr -= x[k]*x[k] - 2*x[k]*y[k] + y[k]*y[k];
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window for one of {variance, std. deviation, z-score}
NumericVector rolling_varsdz(const NumericVector& x, int window, bool pop, CalcTypeVar ctype, bool expanding) {
  
  if (pop && ctype!=ZSCORE) {
    check_x1(x, window);  
  } else {
    check_x(x, window);
  }
  
  int W = window;
  int n_x = x.length();
  int n = 0;
  int var_n = pop ? W : W - 1;
  long double avg = 0, sumsq = 0, delta = 0, varsd = 0, varsd_tmp = 0;
  NumericVector rollx(n_x);
  
  if (expanding) {
    for (int i = 0; i < n_x; ++i) {
      ++n;
      delta = x[i] - avg;
      avg += delta / n;
      sumsq += delta * (x[i] - avg);
      varsd = sumsq / (pop ? n : n - 1);
      varsd_tmp = (ctype == STDDEV || ctype == ZSCORE) ? sqrt(varsd) : varsd;
      if (i >= window - 1) {
        if (ctype == ZSCORE) {
          // handle z-score case to avoid division by zero
          rollx[i] = (i == 0) ? NA_REAL : (x[i] - avg) / varsd_tmp;
        } else {
          rollx[i] = varsd_tmp;  
        }  
      } else {
        rollx[i] = NA_REAL;
      }
    }
    return rollx;
  }
  
  for (int i = 0; i < W; ++i) {
    ++n;
    delta = x[i] - avg;
    avg += delta / n;
    sumsq += delta * (x[i] - avg);
    rollx[i] = NA_REAL;
  }
  
  varsd = sumsq / var_n;
  varsd_tmp = (ctype == STDDEV || ctype == ZSCORE) ? sqrt(varsd) : varsd;
  
  rollx[W - 1] = (ctype == ZSCORE) ? (x[W-1] - avg) / varsd_tmp : varsd_tmp;
  
  long double xi_old = x[0], xi = 0, avg_old = 0;
  
  for (int i = W; i < n_x; ++i) {
    xi = x[i];
    avg_old = avg;
    avg = avg_old + (xi - xi_old) / W;
    varsd += (xi - xi_old)*(xi - avg + xi_old - avg_old) / var_n;
    xi_old = x[i - W + 1];
    
    if (i >= W - 1){
      varsd_tmp = (ctype == STDDEV || ctype == ZSCORE) ? sqrt(varsd) : varsd;
      rollx[i] = (ctype == ZSCORE) ? (xi - avg) / varsd_tmp : varsd_tmp;
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

/*
* PUBLIC EXPORTS
*/

//'
//'Estimates the coefficient (beta) from a univariate linear model y ~ x over a rolling window
//'
//'Estimates \eqn{\hat{\beta_1}} from a ordinary least squares model \eqn{y = \beta_0 + \beta_1 x}.
//'Beta is defined as
//'\deqn{ \code{Beta} = \frac{cov(x, y)}{\sigma^2_x} }
//'where \eqn{cov} is covariance and \eqn{\sigma^2} is variance.
//'
//'@param x A vector or single-column matrix; the independent variable
//'@param y A vector or single-column matrix; the dependent variable
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) covariance, FALSE to calculate sample (n-1) covariance
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingBeta(const NumericVector& x, const NumericVector& y, int window, 
                          bool pop = false, bool expanding = false) {
  return rolling_covcorrbeta(x, y, window, pop, expanding, BETA);
}

//'Calculates simple compounding over a rolling window (e.g. of financial asset return series)
//'
//'Compounding a series of rates of return is defined as
//'\deqn{ \code{Compound} = \big(\prod^n_{i=1} (1+x_i)\big)^{scale / window} - 1 }
//'where \code{scale} is a parameter provided to (optionally) normalize the compounding to a desired period length;
//'this would commonly be used in a financial application where a return could be annualized by (e.g.) setting \code{scale}
//'to 252, 52, 12 or 4 if the time series is daily, weekly, monthly or quarterly, respectively, although any value for \code{scale}
//'may be used.
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param scale Used in the compounding exponent: product^(scale / window) - 1; any number (integer or floating point) may be used
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingCompound(const SEXP& x, int window, long double scale = 1.0, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_compound(R(_, i), window, scale, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates covariance over a rolling window
//'
//'Sample covariance (see parameter \code{pop}) is defined as
//'\deqn{ \code{Cov} = \frac{1}{n-1} \sum^n_{i=1} (x_i - \bar{x})(y_i - \bar{y}) }
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) covariance, FALSE to calculate sample (n-1) covariance
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingCov(const NumericVector& x, const NumericVector& y, int window, 
                         bool pop = false, bool expanding = false) {
  return rolling_covcorrbeta(x, y, window, pop, expanding, COV);
}

//'Calculates Pearson correlation over a rolling window
//'
//'Correlation is defined as
//'\deqn{ \code{Corr} = \frac{cov(x, y)}{\sigma_x\sigma_y} }
//'where \eqn{cov} is the covariance and \eqn{\sigma} is the standard deviation
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) std. deviation, FALSE to calculate sample (n-1) std. deviation
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingCorr(const NumericVector& x, const NumericVector& y, int window, 
                          bool pop = false, bool expanding = false) {
  return rolling_covcorrbeta(x, y, window, pop, expanding, CORR);
}

//'Calculates the excess kurtosis over a rolling window
//'
//'Population excess skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Kurt} = \frac{\frac{1}{n}\sum^n_{i=1}(x_i-\bar{x})^4}{\big(\frac{1}{n}\sum^n_{i=1}(x_i-\bar{x})^2\big)^2} - 3}
//'Sample excess kurtosis (see parameter \code{pop}) is defined as
//'\deqn{ \code{Kurt} = \frac{n(n+1)}{(n-1)(n-2)(n-3)}\frac{\sum^n_{i=1}(x_i-\bar{x})^4}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^2} -3 \frac{(n-1)^2}{(n-2)(n-3)} }
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingKurt(const SEXP& x, int window, bool pop = false, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_skewkurt(R(_, i), window, pop, KURT, expanding);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the maximum over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMax(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_minmax(R(_, i), window, false, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the arithmetic mean over a rolling window
//'
//'The arithmetic mean is defined as
//'\deqn{\code{Mean} = \frac{1}{n}\sum^n_{i=1} x_i}
//'
//'@param x A vector or single-column matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMean(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_mean(R(_, i), window, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the mean absolute error over a rolling window
//'
//'Mean absolute error is defined as
//'\deqn{ \code{MAE} = \frac{1}{n}\sum^n_{i=1}|x_i - y_i| }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingMAE(const NumericVector& x, const NumericVector& y, int window, bool expanding = false) {
  return rolling_meanabserr(x, y, window, expanding);
}

//'Calculates the mean of squares over a rolling window
//'
//'Mean of squares is defined as
//'\deqn{ \code{MS} = \frac{1}{n}\sum^n_{i=1} x^2_i }
//'
//'@param x A vector or single-column matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingMS(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_sumsq(R(_, i), window, expanding, MEAN);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the mean squared error over a rolling window
//'
//'Mean squared error is defined as
//'\deqn{ \code{MSE} = \frac{1}{n}\sum^n_{i=1}(x_i - y_i)^2 }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingMSE(const NumericVector& x, const NumericVector& y, int window, bool expanding = false) {
  return rolling_sqerr(x, y, window, expanding, MEAN);
}

//'Calculates the median over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMedian(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_med(R(_, i), window, expanding);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the minimum over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMin(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for(int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_minmax(R(_, i), window, true, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the product over a rolling window
//'
//'Product is defined as
//'\deqn{ \code{Prod} = \prod^n_{i=1} x_i }
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingProd(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_prod(R(_, i), window, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the root mean squared error over a rolling window
//'
//'Root mean squared error is defined as
//'\deqn{ \code{RMSE} = \sqrt{\frac{1}{n}\sum^n_{i=1}(x_i - y_i)^2 } }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingRMSE(const NumericVector& x, const NumericVector& y, int window, bool expanding = false) {
  return rolling_sqerr(x, y, window, expanding, RMSE);
}

//'Calculates the skewness over a rolling window
//'
//'Population skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Skew} = \sqrt{n}\frac{\sum^n_{i=1}(x_i-\bar{x})^3}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^{3/2}}}
//'Sample skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Skew} = \frac{n\sqrt{n-1}}{n-2}\frac{\sum^n_{i=1}(x_i-\bar{x})^3}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^{3/2}}}
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingSkew(const SEXP& x, int window, bool pop = false, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_skewkurt(R(_, i), window, pop, SKEW, expanding);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the standard deviation over a rolling window
//'
//'Sample standard deviation (see parameter \code{pop}) is defined as
//'\deqn{ \code{Std} = \sqrt{\frac{1}{n-1}\sum^n_{i=1}(x_i - \bar{x})^2} }
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingStd(const SEXP& x, int window, bool pop = false, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_varsdz(R(_, i), window, pop, STDDEV, expanding);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the sum over a rolling window
//'
//'Sum is defined as
//'\deqn{ \code{Sum} = \sum^n_{i=1} x_i }
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingSum(const SEXP& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_sum(R(_, i), window, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the sum product (dot product) over a rolling window
//'
//'Sum product is defined as
//'\deqn{ \code{Sumprod} = \sum^n_{i=1} x_i y_i }
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingSumprod(const NumericVector& x, const NumericVector& y, int window, bool expanding = false) {
  return rolling_sumprod(x, y, window, expanding);
}

//'Calculates the sum of squares over a rolling window
//'
//'Sum of squares is defined as
//'\deqn{ \code{SS} = \sum^n_{i=1} x^2_i }
//'
//'@param x A vector or single-column matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingSS(const NumericVector& x, int window, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_sumsq(R(_, i), window, expanding, SUM);
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the sum of squared errors over a rolling window
//'
//'Sum of squared errors is defined as
//'\deqn{ \code{SSE} = \sum^n_{i=1}(x_i - y_i)^2 }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
//'
//'@param x A vector or single-column matrix
//'@param y A vector or single-column matrix; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericVector RollingSSE(const NumericVector& x, const NumericVector& y, int window, bool expanding = false) {
  return rolling_sqerr(x, y, window, expanding, SUM);
}

//'Calculates the variance over a rolling window
//'
//'Sample variance (see parameter \code{pop}) is defined as
//'\deqn{ \code{Var} = \frac{1}{n-1}\sum^n_{i=1}(x_i - \bar{x})^2 }
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingVar(const SEXP& x, int window, bool pop = false, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_varsdz(R(_, i), window, pop, VAR, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

//'Calculates the z-score (standardized value) over a rolling window
//'
//'z-score is defined as
//'\deqn{ \code{Zscore} = \frac{x - \bar{x}}{\sigma} }
//'where \eqn{\sigma} is the standard deviation of x.
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param pop TRUE to calculate population (n) standard deviation, FALSE to calculate sample (n-1) standard deviation
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingZscore(const SEXP& x, int window, bool pop = false, bool expanding = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix rollx(R.nrow(), R.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    rollx(_, i) = rolling_varsdz(R(_, i), window, pop, ZSCORE, expanding);  
  }
  setnames(R, rollx);
  return rollx;
}

