#ifndef GEANT_RQUEUE
#define GEANT_RQUEUE

#include <deque>
#include "TCondition.h"
#include "TMutex.h"

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_VRUNNABLE
#include "GeantVRunnable.h"
#endif

//______________________________________________________________________________
//   GeantRqueue - A general queue for runnables. 
//______________________________________________________________________________

//______________________________________________________________________________
class GeantRqueue : public TObject
{
private:
   std::deque<GeantVRunnable*> fQueue;      // Queue of runnables
   mutable TMutex              fMutex;      // Mutex for the queue
   TCondition                  fCondition;  // Condition to take next runnable
   Int_t                       fMax;        // Maximum allowed size
   Int_t                       fN;          // Numer of objects in the queue
   Int_t                       fNp;         // Number of prioritized objects

public:

   // Blocking push. The caller is quarantined if the queue is full
   void            Push(GeantVRunnable *data, Bool_t priority=kFALSE);   

   // Try push. If the queue is full return false and the data obkect is not injected
   Bool_t          TryPush(GeantVRunnable *data, Bool_t priority=kFALSE);
   
   // Run next runnable in the queue, possibly. Blocking if the queue is empty.
   GeantVRunnable *RunNext();

   // Try to run next runnable in the queue. Non-blocking even if the queue is empty
   GeantVRunnable *TryRunNext();
   
   // Extract objects from the queue
   
   
