#ifndef GEANT_SYNCOBJECTS
#define GEANT_SYNCOBJECTS
#include <deque>
#include "TCondition.h"
#include "TMutex.h"

using namespace std;

class TStopwatch;
class TObject;

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
   deque<TObject*>   the_queue;
   mutable TMutex    the_mutex;
   TCondition        the_condition_variable;
   TimeCounter      *the_counter;
   Int_t             nobjects;
   Int_t             npriority;
public:
   concurrent_queue(bool counter=false);
   ~concurrent_queue();
   
   Int_t              assigned_workers() const {return the_counter->nthreads;}
   void               push(TObject *data, Bool_t priority=false);
   Int_t              size() const;
   Int_t              size_async() const {return the_queue.size();}
   Bool_t             empty() const;
   Bool_t             empty_async() const {return the_queue.empty();}
   TObject*           wait_and_pop();
   TObject*           wait_and_pop_max(UInt_t nmax, UInt_t &n, TObject **array);
   void               pop_many(UInt_t n, TObject **array);
   Int_t              size_priority() const {return npriority;}
   Int_t              size_objects() const {return nobjects;}
   void               Print();
};
#endif
