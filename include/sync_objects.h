#ifndef GEANT_SYNCOBJECTS
#define GEANT_SYNCOBJECTS
#include <queue>
#include "TCondition.h"
#include "TMutex.h"

class TStopwatch;
class GeantBasket;

//______________________________________________________________________________
struct TimeCounter {
   Int_t       nthreads;
   TStopwatch *timer;
   Double_t    stamp;
   Double_t    realtime[100];
   
   TimeCounter(bool measure_time);
   ~TimeCounter();
   void operator++();
   void operator--();
   void Print();
};         

//______________________________________________________________________________
class concurrent_queue
{
private:
   std::queue<GeantBasket*> the_queue;
   mutable TMutex    the_mutex;
   TCondition        the_condition_variable;
   TimeCounter      *the_counter;
public:
   concurrent_queue(bool counter=false);
   ~concurrent_queue();
   
   Int_t              assigned_workers() const {return the_counter->nthreads;}
   void               push(GeantBasket *data);
   Int_t              size() const {return the_queue.size();}
   Bool_t             empty() const;
   GeantBasket*       try_pop();
   GeantBasket*       wait_and_pop();
   void               Print();
};
#endif
