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

#include <queue>
#include "expanding_median.h"

ExpandingMedian::ExpandingMedian() : hiSize(0), loSize(0) { }

ExpandingMedian::~ExpandingMedian() { }

double ExpandingMedian::GetMedian() {
  if (hiSize > loSize) {
    return hi.top();
  } else if (hiSize < loSize){
    return lo.top();
  } else {
    return (lo.top() + hi.top()) / 2.0;
  }
}

void ExpandingMedian::Insert(double x) {
  if (loSize == 0) {
    lo.push(x);
    ++loSize;
    return;
  }
  if (hiSize == 0) {
    if (x >= lo.top()) {
      hi.push(x);
    } else {
      hi.push(lo.top());
      lo.pop();
      lo.push(x);
    }
    ++hiSize;
    return;
  }
  if (hiSize >= loSize) {
    if (x <= hi.top()) {
      lo.push(x);
    } else {			
      lo.push(hi.top());
      hi.pop();
      hi.push(x);
    }
    ++loSize;
  } else {
    if (x >= lo.top()) {
      hi.push(x);
    } else {			
      hi.push(lo.top());
      lo.pop();
      lo.push(x);
    }
    ++hiSize;
  }	
}
