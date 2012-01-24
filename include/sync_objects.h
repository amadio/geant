#ifndef GEANT_SYNCOBJECTS
#define GEANT_SYNCOBJECTS
#include <queue>
#include "TCondition.h"
#include "TStopwatch.h"
#include "TMutex.h"

class GeantVolumeBasket;

//______________________________________________________________________________
struct TimeCounter {
   Int_t    nthreads;
   TStopwatch timer;
   Double_t stamp;
   Double_t  realtime[100];
   
   TimeCounter();
   void operator++();
   void operator--();
   void Print();
};         

//______________________________________________________________________________
class concurrent_queue
{
private:
    std::queue<GeantVolumeBasket*> the_queue;
    mutable TMutex    the_mutex;
    TCondition        the_condition_variable;
    TimeCounter      *the_counter;
public:
   concurrent_queue(bool counter=false);
   
   void               push(GeantVolumeBasket *data);
   bool               empty() const;
   GeantVolumeBasket* try_pop();
   GeantVolumeBasket* wait_and_pop();
   void               Print();
};
#endif
