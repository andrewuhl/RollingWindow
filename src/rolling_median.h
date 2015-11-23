#ifndef ONLINE_ROLL_MED_H
#define ONLINE_ROLL_MED_H

#include <stdlib.h>

#define inline

typedef double Item;
typedef struct Mediator_t
{
  Item* data;  //circular queue of values
  int*  pos;   //index into `heap` for each value
  int*  heap;  //max/median/min heap holding indexes into `data`.
  int   N;     //allocated size.
  int   idx;   //position in circular queue
  int   minCt; //count of items in min heap
  int   maxCt; //count of items in max heap
} Mediator;

inline int mmless(Mediator* m, int i, int j);
int mmexchange(Mediator* m, int i, int j);
inline int mmCmpExch(Mediator* m, int i, int j);
void minSortDown(Mediator* m, int i);
void maxSortDown(Mediator* m, int i);
inline int minSortUp(Mediator* m, int i);
inline int maxSortUp(Mediator* m, int i);
Mediator* MediatorNew(int nItems);
void MediatorInsert(Mediator* m, Item v);
Item MediatorMedian(Mediator* m);

#endif