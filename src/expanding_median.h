#ifndef EXPANDING_MED_H
#define EXPANDING_MED_H

#include <queue>

class ExpandingMedian {
  
public:	
  double GetMedian();
  void Insert(double);
  ExpandingMedian();
  ~ExpandingMedian();
private:
  int hiSize;
  int loSize;
  std::priority_queue<double> lo;
  std::priority_queue<double, std::vector<double>, std::greater<double> > hi;
};

#endif
