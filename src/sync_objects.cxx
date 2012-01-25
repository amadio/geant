#include "sync_objects.h"
#include "TThread.h"
#include "TGeoManager.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "GeantVolumeBasket.h"

//______________________________________________________________________________
TimeCounter::TimeCounter()
{
// ctor
      timer.Start();
      timer.Stop();
      nthreads = 0;
      stamp = timer.RealTime();
      memset(realtime, 0, 100*sizeof(Double_t));
}

//______________________________________________________________________________
void TimeCounter::operator++() {
// Increment
   TThread::Lock();
   timer.Stop();
   Double_t now = timer.RealTime();
   realtime[nthreads] += now-stamp;
   nthreads++;
   stamp = now;
//      Printf("%d: %f", nthreads, stamp);
   timer.Continue();
   TThread::UnLock();
}   

//______________________________________________________________________________
void TimeCounter::operator--() {
// Decrement
   TThread::Lock();
   timer.Stop();
   Double_t now = timer.RealTime();
   if (nthreads>0) {
      realtime[nthreads] += now-stamp;
      nthreads--;
      stamp = now;
   }   
//      Printf("%d: %f", nthreads, stamp);
   timer.Continue();
   TThread::UnLock();
}

//______________________________________________________________________________
void TimeCounter::Print()
{
// Draw timing statistics.
   timer.Stop();
   Int_t npoints = 0;
   Double_t sum = 0.;
   Int_t i;
   for (i=0; i<100; i++) {
      if (realtime[i]<0.00001) continue;
      npoints++;
      sum += realtime[i];
   }
   Printf("== Percentage of time spent with N threads running (0 = main thread) ==");
   for (i=0; i<npoints; i++) {
      realtime[i] /= sum;
      Printf("%d:  %f%%", i, 100.*realtime[i]);
   }   
   TCanvas *c1 = new TCanvas("c1","Time spent in workers",200,10,700,500);
   c1->SetFillColor(42);
   c1->SetGrid();
   gStyle->SetHistMinimumZero();
   TH1F *h1 = new TH1F("h1","Relative real time spent in workers",npoints,0,npoints);
   h1->SetFillColor(4);
   h1->SetBarWidth(0.4);
   h1->SetBarOffset(0.1);
   h1->SetStats(0);
   h1->SetMinimum(0);
   h1->SetMaximum(1);
   h1->Fill("MT", realtime[0]);
   h1->GetXaxis()->SetBinLabel(1,(new TString("MT"))->Data());
   for (i=2; i<=npoints; i++) {
      TString *s = new TString(Form("%d",i-1));
      h1->Fill(s->Data(), realtime[i-1]);
      h1->GetXaxis()->SetBinLabel(i,s->Data());
   }
   h1->Draw("b");
}

//______________________________________________________________________________
concurrent_queue::concurrent_queue(bool counter)
                 :the_queue(),
                  the_mutex(),
                  the_condition_variable(&the_mutex),
                  the_counter(0)
{
// Concurrent queue constructor.
   if (counter) the_counter = new TimeCounter();
}
                     
//______________________________________________________________________________
void concurrent_queue::push(GeantVolumeBasket *data) {
   the_mutex.Lock();
   the_queue.push(data);
//      Printf("PUSHED basket %s", data->GetName());
   the_condition_variable.Signal();
   the_mutex.UnLock();
}

//______________________________________________________________________________
bool concurrent_queue::empty() const {
   the_mutex.Lock();
   bool is_empty = the_queue.empty();
   the_mutex.UnLock();
   return is_empty;
}

//______________________________________________________________________________
GeantVolumeBasket *concurrent_queue::try_pop() {
   the_mutex.Lock();
   if(the_queue.empty()) {
      the_mutex.UnLock();
      return 0;
   }
   GeantVolumeBasket *popped_value=the_queue.front();
   the_queue.pop();
   the_mutex.UnLock();
   return popped_value;
}

//______________________________________________________________________________
GeantVolumeBasket *concurrent_queue::wait_and_pop() {
   if (the_counter) --(*the_counter);
   the_mutex.Lock();
   while(the_queue.empty()) {
//         Printf("WAITING");
      the_condition_variable.Wait();
   }
        
   GeantVolumeBasket *popped_value=the_queue.front();
   the_queue.pop();
//   Printf("Popped basket %s", popped_value->GetName());
   if (the_counter) ++(*the_counter);
   the_mutex.UnLock();
   return popped_value;
}

//______________________________________________________________________________
void concurrent_queue::Print()
{
   if (the_counter) the_counter->Print();
}   

//______________________________________________________________________________
struct BlockStart {
   Int_t               nthreads;     // number of threads
   Int_t               count;        // counter
   TMutex             *lock;         // Lock for the buffer
   TCondition         *notFull;      // Condition for buffer not full
   TCondition         *notEmpty;     // Condition for buffer not empty

   BlockStart(Int_t num) {
      nthreads = num;
      count    = 0;
      lock     = new TMutex();
      notFull  = new TCondition(lock);
      notEmpty = new TCondition(lock);
   }

   void Start() {
     // Start once then wait N clients to receive.
     // Called only by the main thread
     lock->Lock();
     while (count > 0) {
        Printf("Main thread waiting to give the start...");
        // since already passed once...
        notFull->Wait();
     }
     // Open the door for N threads
     count += nthreads;
     Printf("Main thread waking-up workers");
     notEmpty->Broadcast();
     lock->UnLock();
   }

   void StartN() {
     // Start N times to recieve once
     lock->Lock();
     Int_t tid = TGeoManager::ThreadId();
     // Allow nthreads to pass and increase the counter
     while (count == nthreads) {
        Printf("Thread %d waiting in StartN...", tid);
        notFull->Wait();
     }
     count++;
     notEmpty->Broadcast();
     Printf("(%d) StartN given", tid);
     lock->UnLock();
   }

   void Receive() {
   // Wait until receiving N clients
     lock->Lock();
     Int_t tid = TGeoManager::ThreadId();
     while (count == 0) {
        Printf("Thread %d waiting in Receive...", tid);
        notEmpty->Wait();
     }   
     count--;
     notFull->Broadcast();
     Printf("(%d) Receive.", tid);
     lock->UnLock();
   }

   void ReceiveN() {
   // Wait until receiving N clients
     lock->Lock();
     while (count < nthreads) {
        Printf("Main thread waiting in ReceiveN...");
        notEmpty->Wait();
     }   
     count -= nthreads;
     Printf("ReceiveN broadcasting");
     notFull->Broadcast();
     lock->UnLock();
   }
};

//______________________________________________________________________________
struct ThreadsGate {
   Int_t    ngate;     // Number of threads at the gate
   Int_t    nthreads;  // Total number of threads
   TMutex   mutex;     //
   TMutex   block;     //
   TCondition condvar; // Condition
   TCondition last;    //
   TCondition *condmain;   // Main condition to release on exit
   
   ThreadsGate(Int_t nth, TCondition *cmain=0) : ngate(0), nthreads(nth), mutex(), block(), condvar(&mutex), last(&mutex), condmain(cmain) {}
   ~ThreadsGate() {ngate = nthreads = 0;}
   void              Lock();
   void              UnLock();
   
   void              Sync();
};   

//______________________________________________________________________________
void ThreadsGate::Lock()
{
// Locks the gate.
   block.Lock();
}   

//______________________________________________________________________________
void ThreadsGate::UnLock()
{
// Locks the gate.
   block.UnLock();
}   

//______________________________________________________________________________
void ThreadsGate::Sync()
{
// Enter the gate
   if (nthreads<2) return;        // trivial case
   const char *gateId = (condmain==0)?"start":"stop";
   Int_t tid = TGeoManager::ThreadId();
   block.Lock();                  // lock the block -- new  threads sleep here
   mutex.Lock();                  // lock the mutex
   if (++ngate < nthreads) {      // are we the last thread entering the gate ?
      block.UnLock();             //    not the last one -> unlock the block and
      condvar.Wait();             //    go to sleep 
   } else {                       // yes, we're last
      Printf("-> %s GATE ENTERED BY ALL", gateId);
      condvar.Broadcast();        //    wake everybody up and
      last.Wait();                //    go to sleep till they're all awake... then
      Printf("== %s last exited", gateId);
      if (condmain) {
         Printf("== (%d %s) waking up main thread", tid, gateId);
//         gCurrentBasket->Clear();
         condmain->Broadcast(); // wake up main thread and
      }
      block.UnLock();             //    release the block
   }
   if (--ngate == 1) {            // last but one out ?
      last.Broadcast();           //    yes, wake up last one
   }
   if (!ngate) Printf("%s GATE LEFT BY ALL", gateId);
   mutex.UnLock();
}
