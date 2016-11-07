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
//#include <functional>

#include <Rcpp.h>
#include "rolling_median.h"
#include "expanding_median.h"

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]


/*
* TODO:
* 
*  CACHE LOCALITY
*  ->  any way to improve cache locality?
*  
*  WINDOW SIZE OPTIMIZATION
*  ->  for which window sizes is it optimal to use a naive algorithm vs. single pass?
*  
*  PARALLELIZE
*  ->  add parallel support for:
*        (1) single vector calculations with multi-column matrices, and/or
*        (2) all calculations with vector length > N (N to be determined)
*/

// types of calculations for functions with multiple possible calculations
enum CalcType {COMPOUND, STDDEV, VAR, ZSCORE, BETA, COV, CORR, SKEW, KURT, MEAN, RMSE, SUM, MEANSQ, SUMSQ, MAE, MSE, MEDIAN, SUMPROD, SSE, MIN, MAX, PROD};

// basic error constants
enum ErrorType {XY_LEN_NE, X_LEN_EQ_0, X_LEN_LT_WIN, X_LEN_EQ_1, WIN_LEN_EQ_0, WIN_LEN_EQ_1, WIN_LEN_EQ_N};

// sanity checking function for data input
enum CheckFun {X, X1, X2, X3, X4, XY, XY1};

// stores function arguments for non-data
struct Args {
  int window;
  bool expanding;
  bool pop;
  CalcType ctype;
  std::string na_method;
  double scale;
  CheckFun check_fun;
};

// print errors to R console
void print_error(ErrorType etype, int minwin = 1) {
  std::string errmsg = "";
  switch(etype){
  case XY_LEN_NE:
    errmsg = "vector lengths must be equal";
    break;
  case X_LEN_LT_WIN:
    errmsg = "vector length(s) must be >= window size";
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
    errmsg = "unspecified error";
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
bool check_xn(const NumericVector& x, int window, int minwin){
  if (x.length() < 2){
    print_error(X_LEN_EQ_1);
    return false;
  }
  if (window < minwin){
    print_error(WIN_LEN_EQ_N, minwin);    
    return false;
  }
  if (x.length() < window){
    print_error(X_LEN_LT_WIN);
    return false;
  }
  return true;
}

// check x for cases when length(x) must be >= 2
bool check_x(const NumericVector& x, int window){
  if (x.length() < 2){
    print_error(X_LEN_EQ_1);
    return false;
  }
  if (window < 2){
    print_error(WIN_LEN_EQ_1);
    return false;
  }
  if (x.length() < window){
    print_error(X_LEN_LT_WIN);
    return false;
  }
  return true;
}

// check x for cases when length(x) must be >= 1
bool check_x1(const NumericVector& x, int window){
  if (x.length() == 0){
    print_error(X_LEN_EQ_0);
    return false;
  }
  if (window == 0){
    print_error(WIN_LEN_EQ_0);
    return false;
  }
  if (x.length() < window){
    print_error(X_LEN_LT_WIN);
    return false;
  }
  return true;
}

// check x and y for two-vector calculations
bool check_xy(const NumericVector& x, const NumericVector& y, int window){
  if (x.length() != y.length()){
    print_error(XY_LEN_NE);
    return false;
  }
  if (window < 2){
    print_error(WIN_LEN_EQ_1);
    return false;
  }
  if (x.length() < 2){
    print_error(X_LEN_EQ_1);
    return false;
  }
  if (x.length() < window){
    print_error(X_LEN_LT_WIN);
    return false;
  }
  return true;
}

// check x and y for two-vector calculations
bool check_xy1(const NumericVector& x, const NumericVector& y, int window){
  if (x.length() != y.length()){
    print_error(XY_LEN_NE);
    return false;
  }
  if (window == 0){
    print_error(WIN_LEN_EQ_0);
    return false;
  }
  if (x.length() < window){
    print_error(X_LEN_LT_WIN);
    return false;
  }
  return true;
}

// check x based on check fun type
bool check_function_xy(Args a, const NumericVector&x, const NumericVector&y){
  CheckFun f = a.check_fun;
  if(f == XY){
    return check_xy(x, y, a.window);
  } else if(f == XY1) {
    return check_xy1(x, y, a.window);
  }
  return false;
}

// check x and y based on check fun type
bool check_function_x(Args a, const NumericVector&x){
  CheckFun f = a.check_fun;
  if(a.check_fun == X){
    return check_x(x, a.window);
  } else if(f == X1){
    return check_x1(x, a.window);
  } else if(f == X2){
    return check_xn(x, a.window, 2);
  } else if(f == X3){
    return check_xn(x, a.window, 3);
  } else if(f == X4){
    return check_xn(x, a.window, 4);
  }
  return false;
}

// check if vector is of type NA
bool isna(double x){
  return NumericVector::is_na(x);
}

// substitute missing values for provided replacement value
NumericVector na_sub(const NumericVector& x, double repl, bool last_obs, int maxgap, bool leading) {
  
  if(maxgap < -1){
    stop("maxgap must range between -1 and +Inf");
  }
  
  int n  = x.length();
  NumericVector y(n);
  
  int gap_start_idx = -1, gap_count = 0;
  bool start = false, na_sequence = false;
  double prev_val;
  prev_val = NA_REAL;
  
  if(maxgap == -1){
    maxgap = 2147483647; // largest signed integer
  }
  
  for (int i = 0; i < n; ++i) {
    y[i] = x[i];
    if (NumericVector::is_na(x[i])) {
      if (start | leading) {
        na_sequence = true;
      }
      if (na_sequence) {
        ++gap_count;
        if (gap_start_idx == -1) {
          gap_start_idx = i;
        }
        if (gap_count <= maxgap) {
          if(last_obs){
            if(!(i == gap_count - 1 && leading)){
              repl = prev_val;  
            }
          }
          y[i] = repl;
        } else if (gap_count == maxgap + 1) {
          // we've gone over maxgap, so we need to go back to
          // the start of the gap sequence and replace repl values
          // with NAs
          for (int k = gap_start_idx; k < i; ++k) {
            y[k] = NA_REAL;
          } 
        } // for case gap_count > maxgap + 1, we leave NAs intact.
      }
    } else {
      start = true;
      na_sequence = false;
      gap_start_idx = -1;
      gap_count = 0;
    }
    if(!isna(x[i])){
      prev_val = x[i];  
    }
  }  
  return y;
}

NumericMatrix ok_item_markers(const NumericVector& x, int window, bool expanding) {
  
  // first column holds 0/1 flag; 1 if x_i = NA, 0 otherwise
  // second column holds 0/1 flag; 1 if x_i is calculable , 0 otherwise
  NumericMatrix markers(x.length(), 2);
  
  // number of consecutive non-NA values
  int consec_ok = 0;
  
  for (int i = 0; i < x.length(); ++i) {
    if(isna(x[i])) {
      markers(i, 0) = 1;
      if(!expanding){
        consec_ok = 0;
      }
    } else {
      ++consec_ok;
      if(consec_ok >= window){
        markers(i, 1) = 1;
      }
    }
  }
  return markers;
}

NumericMatrix ok_item_markers_ignore(const NumericVector& x, int window, bool expanding) {
  
  // first column holds 0/1 flag; 1 if x_i = NA, 0 otherwise
  // second column holds 0/1 flag; 1 if x_i is calculable , 0 otherwise
  NumericMatrix markers(x.length(), 2);
  
  for (int i = 0; i < x.length(); ++i) {
    if(isna(x[i])) {
      markers(i, 0) = 1;
    }
    if(i > window - 2){
      markers(i, 1) = 1;  
    }
  }
  return markers;
}

NumericMatrix ok_item_markers_2(const NumericVector& x, const NumericVector& y, int window, bool expanding) {
  
  // first column holds 0/1 flag; 1 if x_i = NA, 0 otherwise
  // second column holds 0/1 flag; 1 if x_i is calculable , 0 otherwise
  NumericMatrix markers(x.length(), 2);
  
  // number of consecutive non-NA values
  int consec_ok = 0;
  
  for (int i = 0; i < x.length(); ++i) {
    if(isna(x[i]) || isna(y[i])) {
      markers(i, 0) = 1;
      if(!expanding){
        consec_ok = 0;
      }
    } else {
      ++consec_ok;
      if(consec_ok >= window){
        markers(i, 1) = 1;
      }
    }
  }
  return markers;
}

NumericMatrix ok_item_markers_2_ignore(const NumericVector& x, const NumericVector& y, int window, bool expanding) {
  
  // first column holds 0/1 flag; 1 if x_i = NA, 0 otherwise
  // second column holds 0/1 flag; 1 if x_i is calculable , 0 otherwise
  NumericMatrix markers(x.length(), 2);
  
  for (int i = 0; i < x.length(); ++i) {
    if(isna(x[i]) || isna(y[i])) {
      markers(i, 0) = 1;
    }
    if(i > window - 2){
      markers(i, 1) = 1;  
    }
  }
  return markers;
}

// maps na-clean data back to original data dimensions by using "marker" matrix
NumericVector map_items(const NumericMatrix& markers, const NumericVector& rollx, int window) {
  NumericVector x(markers.nrow());
  int k = 0;
  for (int i = 0; i < x.length(); ++i) {
    x[i] = NA_REAL;
    // if original item i is not NA
    if(markers(i, 0) == 0){
      x[i] = rollx[k];
      ++k;
      // if original item i is not calculable
      if(markers(i, 1) == 0){
        x[i] = NA_REAL;
      }
    }
  }
  return x;
}

NumericVector nas(int n){
  NumericVector x(n);
  for(int i = 0; i < n; i++) x[i] = NA_REAL;
  return x;
}

void na_handler_method_error(std::string na_method){
  stop("unrecognized value '" + na_method + "' provided for argument `na_method`. `na_method` must be one of ('none', 'ignore', 'window').");
}

// handles NAs for functions with single input vector
NumericVector na_handler(NumericVector (*f)(const NumericVector&, Args),
                         const NumericVector& x, Args a){
  NumericMatrix markers;
  if(a.na_method == "window"){
    markers = ok_item_markers(x, a.window, a.expanding);
  } else if(a.na_method == "ignore") {
    markers = ok_item_markers_ignore(x, a.window, a.expanding);
  } else {
    na_handler_method_error(a.na_method);
  }
  NumericVector x_clean = na_omit(x);
  if(!check_function_x(a, x_clean)){
    return nas(x.length());
  }
  NumericVector rollx = f(x_clean, a);
  return map_items(markers, rollx, a.window);
}

// handles NAs for functions with X & Y input vectors
NumericVector na_handler_2(NumericVector (*f)(const NumericVector&, const NumericVector&, Args), 
                           const NumericVector& x, const NumericVector& y, Args a){
  NumericMatrix xy_markers;
  if(a.na_method == "window"){
    xy_markers = ok_item_markers_2_ignore(x, y, a.window, a.expanding);
  } else if(a.na_method == "ignore"){
    xy_markers = ok_item_markers_2(x, y, a.window, a.expanding);
  } else {
    na_handler_method_error(a.na_method);
  }
  
  LogicalVector xy_ok(x.length());
  int ok_count = 0;
  for(int i = 0; i < x.length(); ++i){
    if(!isna(x[i]) && !isna(y[i])){
      ++ok_count;
      xy_ok[i] = true;
    }
  }
  
  NumericVector x_clean(ok_count);
  NumericVector y_clean(ok_count);
  
  for(int i = 0, k = 0; i < xy_ok.length(); ++i){
    if(xy_ok[i]){
      x_clean[k] = x[i];
      y_clean[k] = y[i];
      ++k;
    }
  }
  
  if(!check_function_xy(a, x_clean, y_clean)){
    return nas(x.length());
  }
  
  NumericVector rollxy = f(x_clean, y_clean, a);
  
  return map_items(xy_markers, rollxy, a.window);
}

/* ------------------------------
* PRIVATE
* 
* Vector-based calculations
* ------------------------------
*/

// calculates rolling compound (useful for financial asset returns)
NumericVector rolling_compound(const NumericVector& x, Args a) {
  
  check_x1(x, a.window);
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double prod = 1;
  for (int i = 0; i < n; ++i) {
    prod *= 1 + x[i];
    if (i >= a.window - 1){
      if (a.expanding) {
        rollx[i] = pow(prod, a.scale / (long double) (i + 1)) - 1;
      } else {
        rollx[i] = pow(prod, a.scale / (long double) a.window) - 1;
        prod /= (1 + x[i - w1]);  
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window for one of {covariance, correlation, beta}
NumericVector rolling_covcorrbeta(const NumericVector& x, const NumericVector& y, Args a) {
  
  if (a.pop) {
    check_xy1(x, y, a.window); 
  } else {
    check_xy(x, y, a.window);  
  }
  
  int W = a.window;
  int n_xy = x.length();
  int w = 0;
  int pop_n = a.pop ? W : W - 1;
  long double avg_x = 0, sumsq_x = 0, delta_x = 0, var_x = 0, sd_x = 0, sum_x = 0;
  long double avg_y = 0, sumsq_y = 0, delta_y = 0, var_y = 0, sd_y = 0, sum_y = 0;
  long double sum_xy = 0;
  long double cov  = 0;
  
  NumericVector rollxy(n_xy);
  
  if (a.expanding) {
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
      
      pop_n = (a.pop ? w : w - 1);
      var_x = sumsq_x / pop_n;
      var_y = sumsq_y / pop_n;
      if(a.ctype == CORR) {
        sd_x = sqrt(var_x);
        sd_y = sqrt(var_y);  
      }
      cov = (sum_xy - sum_x * sum_y / w) / pop_n;
      
      if (i >= a.window - 1) {
        switch (a.ctype) {
        case BETA : rollxy[i] = cov / var_y; break;
        case CORR : rollxy[i] = cov / (sd_x * sd_y); break;
        case COV  : rollxy[i] = cov; break;
        default : break;
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
  if(a.ctype == CORR) {
    sd_x = sqrt(var_x);
    sd_y = sqrt(var_y);  
  }
  cov = (sum_xy - sum_x * sum_y / W) / pop_n;
  
  switch (a.ctype) {
  case BETA : rollxy[W - 1] = cov / var_y; break;
  case CORR : rollxy[W - 1] = cov / (sd_x * sd_y); break;
  case COV  : rollxy[W - 1] = cov; break;
  default: break;
  }
  
  // std dev terms
  long double xi_old = x[0], xi = 0, avg_old_x = 0;
  long double yi_old = y[0], yi = 0, avg_old_y = 0;
  
  for (int i = W; i < n_xy; ++i) {
    
    // std dev of x
    xi = x[i];
    if (a.ctype == BETA || a.ctype == CORR) {
      avg_old_x = avg_x;
      avg_x = avg_old_x + (xi - xi_old) / W;
      var_x += (xi - xi_old)*(xi - avg_x + xi_old - avg_old_x) / pop_n;
      if(a.ctype == CORR)
        sd_x = sqrt(var_x);
    }
    
    // std dev of y
    yi = y[i];
    if (a.ctype == BETA || a.ctype == CORR) {
      avg_old_y = avg_y;
      avg_y = avg_old_y + (yi - yi_old) / W;
      var_y += (yi - yi_old)*(yi - avg_y + yi_old - avg_old_y) / pop_n;
      if(a.ctype == CORR)
        sd_y = sqrt(var_y);
    }
    
    // cov of x,y
    sum_xy += xi * yi - xi_old * yi_old;
    sum_x += xi - xi_old;
    sum_y += yi - yi_old;
    cov = (sum_xy - sum_x * sum_y / W) / pop_n;
    
    switch (a.ctype) {
    case BETA : rollxy[i] = cov / var_y; break;
    case CORR : rollxy[i] = cov / (sd_x * sd_y); break;
    case COV  : rollxy[i] = cov; break;
    default: break;
    }
    
    xi_old = x[i - W + 1];
    yi_old = y[i - W + 1];
  }
  return rollxy;
}

// calculates rolling mean
NumericVector rolling_mean(const NumericVector& x, Args a) {
  
  check_x1(x, a.window);
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += x[i];
    if (i >= a.window - 1){
      if (a.expanding) {
        rollx[i] = sum / (i + 1);
      } else {
        rollx[i] = sum / a.window;
        sum -= x[i - w1];  
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling mean absolute error
NumericVector rolling_meanabserr(const NumericVector& x, const NumericVector& y, Args a) {
  
  check_xy1(x, y, a.window);
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sumabserr = 0;
  for (int i = 0; i < n; ++i) {
    sumabserr += std::abs(x[i] - y[i]);
    if (i >= a.window - 1) {
      if (a.expanding) {
        rollx[i] = sumabserr / (i + 1);
      } else {
        rollx[i] = sumabserr / a.window;  
        sumabserr -= std::abs(x[i-w1] - y[i - w1]);
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calcuates rolling median
NumericVector rolling_med(const NumericVector& x, Args a) {
  
  check_x1(x, a.window);
  
  int n = x.length();
  NumericVector rollx(n);
  
  if (a.expanding) {
    
    ExpandingMedian* median = new ExpandingMedian();
    
    for (int i = 0; i < n; ++i) {
      median->Insert(x[i]);
      if (i >= a.window - 1) {
        rollx[i] = median->GetMedian();  
      } else {
        rollx[i] = NA_REAL;
      }
    }
    delete median;
    
  } else {
    
    Mediator* m = MediatorNew(a.window);
    
    for (int i = 0; i < n; ++i) {
      MediatorInsert(m, x[i]);
      if (i >= a.window - 1) {
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
NumericVector rolling_minmax(const NumericVector& x, Args a) {
  
  check_x1(x, a.window);
  
  int n  = x.length();
  NumericVector rollx(n);
  
  if (a.expanding) {
    long double minmax = x[0];
    for (int i = 0; i < n; ++i) {
      if (a.ctype == MIN) {
        if (x[i] < minmax) {
          minmax = x[i];
        }
      } else {
        if (x[i] > minmax) {
          minmax = x[i];
        }
      }
      if (i >= a.window - 1) {
        rollx[i] = minmax;  
      } else {
        rollx[i] = NA_REAL;
      }
    }
    return rollx;
  }
  
  std::deque< std::pair<long double, int> > deck;
  for (int i = 0; i < x.size(); ++i) {
    if(a.ctype == MIN) {
      while (!deck.empty() && deck.back().first >= x[i])
        deck.pop_back();
    } else {
      while (!deck.empty() && deck.back().first <= x[i])
        deck.pop_back();  
    }
    deck.push_back(std::make_pair(x[i], i));
    
    while(deck.front().second <= i - a.window)
      deck.pop_front();
    
    long double min = deck.front().first;
    if (i < a.window - 1) {
      rollx[i] = NA_REAL;
    } else {
      rollx[i] = min;
    }     
  }
  return rollx;
}

// calculates rolling product
NumericVector rolling_prod(const NumericVector& x, Args a) {
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double prod = 1;
  for (int i = 0; i < n; ++i) {
    prod *= x[i];
    if (i >= a.window - 1) {
      rollx[i] = prod;
      if (!a.expanding)
        prod /= x[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window for one of {skewness, kurtosis}
NumericVector rolling_skewkurt(const NumericVector& x, Args a) {
  
  // ----------------------------
  // expanding window calculation  
  // ----------------------------
  
  int W = a.window;
  int n_x = x.length();
  NumericVector rollx(n_x);
  
  if (a.expanding) {
    
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
      if (a.ctype == KURT) {
        M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;  
      }
      M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
      M2 += term1;
      // ----------- END ATTRIBUTION ----------------------------
      
      if (i < W - 1) {
        rollx[i] = NA_REAL;
      } else {
        if (a.ctype == SKEW) {
          pop_skew = a.pop ? double(n)*sqrt(double(n)-1)/(double(n)-2): sqrt(double(n));
          rollx[i] = pop_skew * M3 / pow(M2, 1.5);
        } else if (a.ctype == KURT) {
          rollx[i] = n * M4 / (M2*M2) - 3;
          if(!a.pop){
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
    
    if (a.ctype == SKEW) {
      pop_skew = a.pop ? double(W)*sqrt(double(W)-1)/(double(W)-2): sqrt(double(W));
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
        
        if (a.ctype == SKEW) {
          M3 = sumx3 - c1_M3*sumx*sumx2 + c2_M3*pow(sumx, 3);
          rollx[i] = pop_skew * M3 / pow(M2, 1.5);  
        } else if (a.ctype == KURT) {        
          M4 = sumx4 + sumx*(c1_M4*sumx*sumx2 + c2_M4*sumx3 + c3_M4*pow(sumx, 3));
          rollx[i] = W*M4 / (M2*M2) - 3;
          if (!a.pop) {
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
NumericVector rolling_sum(const NumericVector& x, Args a) {
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += x[i];
    if (i >= a.window - 1) {
      rollx[i] = sum;
      if (!a.expanding)
        sum -= x[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling sum product
NumericVector rolling_sumprod(const NumericVector& x, const NumericVector& y, Args a) {
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sumprod = 0;
  for (int i = 0; i < n; ++i) {
    sumprod += x[i] * y[i];
    if (i >= a.window - 1) {
      rollx[i] = sumprod;
      if (!a.expanding)
        sumprod -= x[i - w1] * y[i - w1];
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates one of rolling {mean of squares, sum of squares}
NumericVector rolling_sumsq(const NumericVector& x, Args a) {
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sumsq = 0;
  for (int i = 0; i < n; ++i) {
    sumsq += x[i] * x[i];
    if (i >= a.window - 1) {
      rollx[i] = sumsq; 
      if (a.ctype == MEANSQ) {
        rollx[i] /= a.expanding ? i + 1 : a.window;
      }
      if (!a.expanding) {
        sumsq -= x[i - w1] * x[i - w1];
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling sum of squared errors
NumericVector rolling_sqerr(const NumericVector& x, const NumericVector& y, Args a) {
  
  int n  = x.length(), w1 = a.window - 1;
  NumericVector rollx(n);
  long double sumsqerr = 0;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    sumsqerr += x[i]*x[i] - 2*x[i]*y[i] + y[i]*y[i];
    if (i >= a.window - 1) {
      rollx[i] = sumsqerr;
      if (a.ctype == MSE || a.ctype == RMSE) {
        rollx[i] /= a.expanding ? i + 1 : a.window;
        if (a.ctype == RMSE) {
          rollx[i] = sqrt(rollx[i]);
        }
      }
      if (!a.expanding) {
        k = i - w1;
        sumsqerr -= x[k]*x[k] - 2*x[k]*y[k] + y[k]*y[k];
      }
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// calculates rolling window Sharpe Ratio
NumericVector rolling_sharpe(const NumericVector& x, const NumericVector& Rf, Args a) {
  
  int W = a.window;
  int n_x = x.length();
  int n = 0;
  int var_n = a.pop ? W : W - 1, w1 = a.window - 1;
  long double avg = 0, sumsq = 0, delta = 0, varsd = 0, prod = 1.0;
  NumericVector rollx(n_x);
  
  if (a.expanding) {
    for (int i = 0; i < n_x; ++i) {
      prod *= 1 + x[i] - Rf[i];
      ++n;
      delta = x[i] - avg;
      avg += delta / n;
      sumsq += delta * (x[i] - avg);
      varsd = sqrt(a.scale * sumsq / (n - 1));
      if (i >= a.window - 1 && i > 0) {
        rollx[i] = (pow(prod, a.scale / (long double) (i + 1)) - 1) / varsd;
      } else {
        rollx[i] = NA_REAL;
      }
    }
    return rollx;
  }
  
  for (int i = 0; i < W; ++i) {
    prod *= 1 + x[i] - Rf[i];
    ++n;
    delta = x[i] - avg;
    avg += delta / n;
    sumsq += delta * (x[i] - avg);
    rollx[i] = NA_REAL;
  }
  
  varsd = sumsq / var_n;
  
  rollx[W - 1] = (pow(prod, a.scale / (long double) a.window) - 1) / sqrt(varsd * a.scale);
  prod /= (1 + x[0] - Rf[0]);
  
  long double xi_old = x[0], xi = 0, avg_old = 0;
  
  for (int i = W; i < n_x; ++i) {
    prod *= 1 + x[i] - Rf[i];
    xi = x[i];
    avg_old = avg;
    avg = avg_old + (xi - xi_old) / W;
    varsd += (xi - xi_old)*(xi - avg + xi_old - avg_old) / var_n;
    xi_old = x[i - W + 1];
    if (i >= W - 1){
      rollx[i] = ((pow(prod, a.scale / (long double) a.window) - 1)) / sqrt(varsd * a.scale);
    } else {
      rollx[i] = NA_REAL;  
    }
    prod /= (1 + x[i - w1] - Rf[i - w1]);
  }
  return rollx;
}


// // calculates rolling window Sharpe Ratio
// NumericVector rolling_sharpe(const NumericVector& x, const NumericVector& Rf, Args a) {
//   
//   int W = a.window;
//   int n_x = x.length();
//   int n = 0;
//   int var_n = a.pop ? W : W - 1, w1 = a.window - 1;
//   long double avg = 0, sumsq = 0, delta = 0, varsd = 0, prod = 1.0;
//   NumericVector rollx(n_x);
//   
//   // this allows for the caller to supply Rf as a vector of length 1.
//   // for readability, we create a NumericVector by repeating Rf[0] x.length() times
//   // so that we can later subtract Rf_vec[i] from x[i]
//   NumericVector Rf_vec;
//   if(Rf.length() == 1){
//     Rf_vec = rep(Rf[0], x.length());
//   } else {
//     Rf_vec = Rf;
//   }
//   
//   // for(int i = 0; i < x.length(); i++){
//   //   Rcout << i + 1 << "\t" << x[i] << "\t " << Rf_vec[i] << std::endl;
//   // }
//   
//   if (a.expanding) {
//     for (int i = 0; i < n_x; ++i) {
//       prod *= 1 + x[i] - Rf_vec[i];
//       ++n;
//       delta = x[i] - avg;
//       avg += delta / n;
//       sumsq += delta * (x[i] - avg);
//       varsd = sqrt(a.scale * sumsq / (n - 1));
//       if (i >= a.window - 1 && i > 0) {
//         rollx[i] = (pow(prod, a.scale / (long double) (i + 1)) - 1) / varsd;
//       } else {
//         rollx[i] = NA_REAL;
//       }
//     }
//     return rollx;
//   }
//   
//   for (int i = 0; i < W; ++i) {
//     prod *= 1 + x[i] - Rf_vec[i];
//     ++n;
//     delta = x[i] - avg;
//     avg += delta / n;
//     sumsq += delta * (x[i] - avg);
//     rollx[i] = NA_REAL;
//   }
//   
//   varsd = sumsq / var_n;
//   
//   rollx[W - 1] = (pow(prod, a.scale / (long double) a.window) - 1) / sqrt(varsd * a.scale);
//   prod /= (1 + x[0] - Rf_vec[0]);
//   
//   long double xi_old = x[0], xi = 0, avg_old = 0;
//   
//   for (int i = W; i < n_x; ++i) {
//     prod *= 1 + x[i] - Rf_vec[i];
//     xi = x[i];
//     avg_old = avg;
//     avg = avg_old + (xi - xi_old) / W;
//     varsd += (xi - xi_old)*(xi - avg + xi_old - avg_old) / var_n;
//     xi_old = x[i - W + 1];
//     if (i >= W - 1){
//       rollx[i] = ((pow(prod, a.scale / (long double) a.window) - 1)) / sqrt(varsd * a.scale);
//     } else {
//       rollx[i] = NA_REAL;  
//     }
//     prod /= (1 + x[i - w1] - Rf_vec[i - w1]);
//   }
//   return rollx;
// }

// calculates rolling window for one of {variance, std. deviation, z-score}
NumericVector rolling_varsdz(const NumericVector& x, Args a) {
  
  int W = a.window;
  int n_x = x.length();
  int n = 0;
  int var_n = a.pop ? W : W - 1;
  long double avg = 0, sumsq = 0, delta = 0, varsd = 0, varsd_tmp = 0;
  NumericVector rollx(n_x);
  
  if (a.expanding) {
    for (int i = 0; i < n_x; ++i) {
      ++n;
      delta = x[i] - avg;
      avg += delta / n;
      sumsq += delta * (x[i] - avg);
      varsd = sumsq / (a.pop ? n : n - 1);
      varsd_tmp = (a.ctype == STDDEV || a.ctype == ZSCORE) ? sqrt(varsd) : varsd;
      if (i >= a.window - 1) {
        if (a.ctype == ZSCORE) {
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
  varsd_tmp = (a.ctype == STDDEV || a.ctype == ZSCORE) ? sqrt(varsd) : varsd;
  
  rollx[W - 1] = (a.ctype == ZSCORE) ? (x[W-1] - avg) / varsd_tmp : varsd_tmp;
  
  long double xi_old = x[0], xi = 0, avg_old = 0;
  
  for (int i = W; i < n_x; ++i) {
    xi = x[i];
    avg_old = avg;
    avg = avg_old + (xi - xi_old) / W;
    varsd += (xi - xi_old)*(xi - avg + xi_old - avg_old) / var_n;
    xi_old = x[i - W + 1];
    
    if (i >= W - 1){
      varsd_tmp = (a.ctype == STDDEV || a.ctype == ZSCORE) ? sqrt(varsd) : varsd;
      rollx[i] = (a.ctype == ZSCORE) ? (xi - avg) / varsd_tmp : varsd_tmp;
    } else {
      rollx[i] = NA_REAL;  
    }
  }
  return rollx;
}

// columnwise apply() for function with most basic signature (SEXP, int, bool, bool)
NumericMatrix apply_rolling_FUN_1(NumericVector (*f)(const NumericVector&, Args), 
                                  const SEXP& x, Args a){
  NumericMatrix xmat = as_matrix(x);
  if(xmat.ncol() > 0) {
    if(!check_function_x(a, xmat(_, 0))){
      Rcpp::stop("Not enough data points to calculate the rolling statistic.");
    }
  }
  NumericMatrix rollx(xmat.nrow(), xmat.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    if(a.na_method == "none"){
      rollx(_, i) = f(xmat(_, i), a);
    } else {
      rollx(_, i) = na_handler(f, xmat(_, i), a);
    }
  }
  setnames(xmat, rollx);
  return rollx;
}

// columnwise apply() for function with most basic signature (SEXP, int, bool, bool, CalcType)
NumericMatrix apply_rolling_FUN_2(NumericVector (*f)(const NumericVector&, const NumericVector&, Args), 
                                  const SEXP& x, const SEXP& y, Args a){
  NumericMatrix xmat = as_matrix(x);
  NumericMatrix ymat = as_matrix(y);
  if(ymat.ncol() > 1){
    Rcpp::stop("y must be a vector or contain exactly 1 column of data.");
  }
  NumericVector yvec = ymat(_, 0);
  if(xmat.ncol() > 0) {
    if(!check_function_xy(a, xmat(_, 0), yvec)){
      Rcpp::stop("Not enough data points to calculate the rolling statistic.");
    }
  }
  NumericMatrix rollx(xmat.nrow(), xmat.ncol());
  for (int i = 0; i < rollx.ncol(); ++i) {
    if(a.na_method == "none"){
      rollx(_, i) = f(xmat(_, i), yvec, a);
    } else {
      rollx(_, i) = na_handler_2(f, xmat(_, i), yvec, a);
    }
  }
  setnames(xmat, rollx);
  return rollx;
}


/*
* PUBLIC API
*/

//'Replaces NAs with either the last non-NA observation, or a specified replacement value
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param repl The value used to replace NAs
//'@param last_obs TRUE if the last non-NA observation should be used instead of a specified replacement value
//'@param maxgap Runs of more than \code{maxgap} are retained. By default, the argument is set to 0, which 
//'will preserve all NAs. If set to -1, the argument is set to infinity such that sequences of any size of NAs
//'will be replaced with the last non-NA observation.
//'@param leading If TRUE, leading NAs will be replaced by the value specified by \code{repl}, otherwise leading
//'NAs will be preserved
// [[Rcpp::export]]
NumericMatrix NaSub(const SEXP& x, double repl = NA_REAL, bool last_obs = true, int maxgap = -1, bool leading = false) {
  NumericMatrix R = as_matrix(x);
  NumericMatrix y(R.nrow(), R.ncol());
  for (int i = 0; i < y.ncol(); ++i) {
    y(_, i) = na_sub(R(_, i), repl, last_obs, maxgap, leading);
  }
  setnames(R, y);
  return y;
}

//'Estimates the coefficient (beta) from a univariate linear model y ~ x over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix; the independent variable
//'@param y A vector, matrix, list or other object coercible to a vector; the dependent variable
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic//'
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) covariance, FALSE to calculate sample (n-1) covariance
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Estimates \eqn{\hat{\beta_1}} from an ordinary least squares model \eqn{y = \beta_0 + \beta_1 x}.
//'Beta is defined as
//'\deqn{ \code{Beta} = \frac{cov(x, y)}{\sigma^2_x} }
//'where \eqn{cov} is covariance and \eqn{\sigma^2} is variance.
//'If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column against the vector \code{y}. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericVector RollingBeta(const SEXP& x, const SEXP& y, int window, 
                          bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.pop = pop; a.ctype = BETA; a.check_fun = a.pop ? XY1 : XY;
  return apply_rolling_FUN_2(rolling_covcorrbeta, x, y, a);
}

//'Calculates simple compounding over a rolling window (e.g. of financial asset return series)
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param scale Used in the compounding exponent: product^(scale / window) - 1; any number (integer or floating point) may be used
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Compounding a series of rates of return is defined as
//'\deqn{ \code{Compound} = \big(\prod^n_{i=1} (1+x_i)\big)^{scale / window} - 1 }
//'where \code{scale} is a parameter provided to (optionally) normalize the compounding to a desired period length;
//'this would commonly be used in a financial application where a return could be annualized by (e.g.) setting \code{scale}
//'to 252, 52, 12 or 4 if the time series is daily, weekly, monthly or quarterly, respectively, although any value for \code{scale}
//'may be used.
// [[Rcpp::export]]
NumericMatrix RollingCompound(const SEXP& x, int window, long double scale = 1.0, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = COMPOUND; a.scale = scale; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_compound, x, a);
}

//'Calculates covariance over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) covariance, FALSE to calculate sample (n-1) covariance
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Sample covariance (see parameter \code{pop}) is defined as
//'\deqn{ \code{Cov} = \frac{1}{n-1} \sum^n_{i=1} (x_i - \bar{x})(y_i - \bar{y}) }
// [[Rcpp::export]]
NumericVector RollingCov(const SEXP& x, const SEXP& y, int window, 
                         bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.pop = pop, a.ctype = COV; a.check_fun = a.pop ? XY1 : XY;
  return apply_rolling_FUN_2(rolling_covcorrbeta, x, y, a);
}

//'Calculates Pearson correlation over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) std. deviation, FALSE to calculate sample (n-1) std. deviation
//'determines the starting number of observations used to calculate the statistic
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Correlation is defined as
//'\deqn{ \code{Corr} = \frac{cov(x, y)}{\sigma_x\sigma_y} }
//'where \eqn{cov} is the covariance and \eqn{\sigma} is the standard deviation
// [[Rcpp::export]]
NumericVector RollingCorr(const SEXP& x, const SEXP& y, int window, 
                          bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.pop = pop, a.ctype = CORR; a.check_fun = a.pop ? XY1 : XY;
  return apply_rolling_FUN_2(rolling_covcorrbeta, x, y, a);
}

//'Calculates the excess kurtosis over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Population excess skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Kurt} = \frac{\frac{1}{n}\sum^n_{i=1}(x_i-\bar{x})^4}{\big(\frac{1}{n}\sum^n_{i=1}(x_i-\bar{x})^2\big)^2} - 3}
//'Sample excess kurtosis (see parameter \code{pop}) is defined as
//'\deqn{ \code{Kurt} = \frac{n(n+1)}{(n-1)(n-2)(n-3)}\frac{\sum^n_{i=1}(x_i-\bar{x})^4}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^2} -3 \frac{(n-1)^2}{(n-2)(n-3)} }
// [[Rcpp::export]]
NumericMatrix RollingKurt(const SEXP& x, int window, bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.pop = pop, a.ctype = KURT; a.check_fun = X4;
  return apply_rolling_FUN_1(rolling_skewkurt, x, a);
}

//'Calculates the maximum over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
// [[Rcpp::export]]
NumericMatrix RollingMax(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MAX; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_minmax, x, a);
}

//'Calculates the arithmetic mean over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details The arithmetic mean is defined as
//'\deqn{\code{Mean} = \frac{1}{n}\sum^n_{i=1} x_i}
// [[Rcpp::export]]
NumericMatrix RollingMean(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MEAN; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_mean, x, a);
}

//'Calculates the mean absolute error over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Mean absolute error is defined as
//'\deqn{ \code{MAE} = \frac{1}{n}\sum^n_{i=1}|x_i - y_i| }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
// [[Rcpp::export]]
NumericVector RollingMAE(const SEXP& x, const SEXP& y, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MAE; a.check_fun = XY1;
  return apply_rolling_FUN_2(rolling_meanabserr, x, y, a);
}

//'Calculates the mean of squares over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Mean of squares is defined as
//'\deqn{ \code{MS} = \frac{1}{n}\sum^n_{i=1} x^2_i }
// [[Rcpp::export]]
NumericVector RollingMS(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MEANSQ; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_sumsq, x, a);
}

//'Calculates the mean squared error over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector or single-column matrix; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Mean squared error is defined as
//'\deqn{ \code{MSE} = \frac{1}{n}\sum^n_{i=1}(x_i - y_i)^2 }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
// [[Rcpp::export]]
NumericVector RollingMSE(const SEXP& x, const SEXP& y, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MSE; a.check_fun = XY1;
  return apply_rolling_FUN_2(rolling_sqerr, x, y, a);
}

//'Calculates the median over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMedian(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MEDIAN; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_med, x, a);
}

//'Calculates the minimum over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details If input \code{x} contains more than one column, the rolling window calculation will be performed on each
//'column. The dimensions of the return value will be the same as input \code{x}.
// [[Rcpp::export]]
NumericMatrix RollingMin(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = MIN; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_minmax, x, a);
}

//'Calculates the product over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Product is defined as
//'\deqn{ \code{Prod} = \prod^n_{i=1} x_i }
// [[Rcpp::export]]
NumericMatrix RollingProd(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = PROD; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_prod, x, a);
}

//'Calculates the root mean squared error over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Root mean squared error is defined as
//'\deqn{ \code{RMSE} = \sqrt{\frac{1}{n}\sum^n_{i=1}(x_i - y_i)^2 } }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
// [[Rcpp::export]]
NumericVector RollingRMSE(const SEXP& x, const SEXP& y, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = RMSE; a.check_fun = XY1;
  return apply_rolling_FUN_2(rolling_sqerr, x, y, a);
}

//'Calculates the skewness over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Population skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Skew} = \sqrt{n}\frac{\sum^n_{i=1}(x_i-\bar{x})^3}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^{3/2}}}
//'Sample skewness (see parameter \code{pop}) is defined as
//'\deqn{ \code{Skew} = \frac{n\sqrt{n-1}}{n-2}\frac{\sum^n_{i=1}(x_i-\bar{x})^3}{\big(\sum^n_{i=1}(x_i-\bar{x})^2\big)^{3/2}}}
// [[Rcpp::export]]
NumericMatrix RollingSkew(const SEXP& x, int window, bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = SKEW; a.pop = pop; a.check_fun = a.pop ? X3 : X2;
  return apply_rolling_FUN_1(rolling_skewkurt, x, a);
}

//'Calculates the standard deviation over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Sample standard deviation (see parameter \code{pop}) is defined as
//'\deqn{ \code{Std} = \sqrt{\frac{1}{n-1}\sum^n_{i=1}(x_i - \bar{x})^2} }
// [[Rcpp::export]]
NumericMatrix RollingStd(const SEXP& x, int window, bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = STDDEV; a.pop = pop; a.check_fun = a.pop ? X1 : X;
  return apply_rolling_FUN_1(rolling_varsdz, x, a);
}

//'Calculates the sum over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Sum is defined as
//'\deqn{ \code{Sum} = \sum^n_{i=1} x_i }
// [[Rcpp::export]]
NumericMatrix RollingSum(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = SUM; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_sum, x, a);
}

//'Calculates the sum product (dot product) over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Sum product is defined as
//'\deqn{ \code{Sumprod} = \sum^n_{i=1} x_i y_i }
// [[Rcpp::export]]
NumericVector RollingSumprod(const SEXP& x, const SEXP& y, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = SUMPROD; a.check_fun = XY1;
  return apply_rolling_FUN_2(rolling_sumprod, x, y, a);
}

//'Calculates the sum of squares over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Sum of squares is defined as
//'\deqn{ \code{SS} = \sum^n_{i=1} x^2_i }
// [[Rcpp::export]]
NumericVector RollingSS(const SEXP& x, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = SUMSQ; a.check_fun = X1;
  return apply_rolling_FUN_1(rolling_sumsq, x, a);
}

//'Calculates the sum of squared errors over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param y A vector, matrix, list or other object coercible to a vector; the error to be subtracted from x
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@return A vector; rows 1 to (length(window) - 1) contain NAs.
//'@details Sum of squared errors is defined as
//'\deqn{ \code{SSE} = \sum^n_{i=1}(x_i - y_i)^2 }
//'where \eqn{x} and \eqn{y} represent the predicted and actual (true) values.
// [[Rcpp::export]]
NumericMatrix RollingSSE(const SEXP& x, const SEXP& y, int window, bool expanding = false, std::string na_method = "none") {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = SSE; a.check_fun = XY1;
  return apply_rolling_FUN_2(rolling_sqerr, x, y, a);
}

//'Calculates the variance over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) statistic, FALSE to calculate sample (n-1) statistic
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details Sample variance (see parameter \code{pop}) is defined as
//'\deqn{ \code{Var} = \frac{1}{n-1}\sum^n_{i=1}(x_i - \bar{x})^2 }
// [[Rcpp::export]]
NumericMatrix RollingVar(const SEXP& x, int window, bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = VAR; a.pop = pop; a.check_fun = a.pop ? X1 : X;
  return apply_rolling_FUN_1(rolling_varsdz, x, a);
}

//'Calculates the z-score (standardized value) over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param pop TRUE to calculate population (n) standard deviation, FALSE to calculate sample (n-1) standard deviation
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details z-score is defined as
//'\deqn{ \code{Zscore} = \frac{x - \bar{x}}{\sigma} }
//'where \eqn{\sigma} is the standard deviation of x.
// [[Rcpp::export]]
NumericMatrix RollingZscore(const SEXP& x, int window, bool expanding = false, std::string na_method = "none", bool pop = false) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = ZSCORE; a.pop = pop; a.check_fun = X;
  return apply_rolling_FUN_1(rolling_varsdz, x, a);
}

//'Calculates the Sharpe Ratio over a rolling window
//'
//'@param x A vector, matrix, list or other object coercible to a matrix
//'@param Rf The risk-free rate, which is must be a vector, matrix, list or other object coercible to a vector.
//'Rf must have the same number of observations as x.
//'@param window The size of the rolling window
//'@param expanding TRUE to calculate an expanding window instead of a rolling window; if TRUE, window
//'determines the starting number of observations used to calculate the statistic
//'@param na_method One of "none" (default), "ignore" or "window" (see notes on NA handling)
//'@param scale Used in the compounding exponent: product^(scale / window) - 1; any number (integer or floating point) may be used
//'@return A matrix; rows 1 to (length(window) - 1) contain NAs.
//'@details The Sharpe Ratio is defined as
//'\deqn{ \code{Zscore} = \frac{x - Rf}{\sigma} }
//'where \eqn{\sigma} is the standard deviation of x.
// [[Rcpp::export]]
NumericMatrix RollingSharpe(const SEXP& x, const SEXP& Rf, int window, bool expanding = false, std::string na_method = "none", double scale = 12) {
  Args a; a.window = window; a.expanding = expanding; a.na_method = na_method; a.ctype = ZSCORE; a.scale = scale, a.check_fun = XY;
  return apply_rolling_FUN_2(rolling_sharpe, x, Rf, a);
}

