#ifndef GEANTV_MIC
#include "sync_objects.h"
#include "TThread.h"
#include "TGeoManager.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TStopwatch.h"

//______________________________________________________________________________
TimeCounter::TimeCounter(bool measure_time) : nthreads(0), timer(0), stamp(0) {
  // ctor
  if (measure_time) {
    timer = new TStopwatch();
    timer->Start();
    timer->Stop();
    stamp = timer->RealTime();
    memset(realtime, 0, 100 * sizeof(double));
  }
}

//______________________________________________________________________________
TimeCounter::~TimeCounter() {
  // destructor.
  delete timer;
}

//______________________________________________________________________________
TimeCounter &TimeCounter::operator++() {
  // Increment
  TThread::Lock();
  if (timer) {
    timer->Stop();
    double now = timer->RealTime();
    realtime[nthreads] += now - stamp;
    stamp = now;
    timer->Continue();
  }
  nthreads++;
  //      Printf("%d: %f", nthreads, stamp);
  TThread::UnLock();
  return *this;
}

//______________________________________________________________________________
TimeCounter &TimeCounter::operator--() {
  // Decrement
  TThread::Lock();
  if (timer) {
    timer->Stop();
    double now = timer->RealTime();
    if (nthreads > 0) {
      realtime[nthreads] += now - stamp;
      nthreads--;
      stamp = now;
    }
    timer->Continue();
  } else {
    if (nthreads > 0)
      nthreads--;
  }
  //      Printf("%d: %f", nthreads, stamp);
  TThread::UnLock();
  return *this;
}

//______________________________________________________________________________

void TimeCounter::Print() {
  // Draw timing statistics.
  if (!timer)
    return;
  timer->Stop();
  int npoints = 0;
  double sum = 0.;
  int i;
  for (i = 0; i < 100; i++) {
    if (realtime[i] < 0.00001)
      continue;
    npoints++;
    sum += realtime[i];
  }
  Printf("== Percentage of time spent with N threads running (0 = main thread) ==");
  for (i = 0; i < npoints; i++) {
    realtime[i] /= sum;
    Printf("%d:  %f%%", i, 100. * realtime[i]);
  }

  TCanvas *c1 = new TCanvas("c1", "Time spent in workers", 200, 10, 700, 500);
  c1->SetFillColor(42);
  c1->SetGrid();
  gStyle->SetHistMinimumZero();
  TH1F *h1 = new TH1F("h1", "Relative real time spent in workers", npoints, 0, npoints);
  h1->SetFillColor(4);
  h1->SetBarWidth(0.4);
  h1->SetBarOffset(0.1);
  h1->SetStats(0);
  h1->SetMinimum(0);
  h1->SetMaximum(1);
  h1->Fill("MT", realtime[0]);
  h1->GetXaxis()->SetBinLabel(1, (new TString("MT"))->Data());
  for (i = 2; i <= npoints; i++) {
    TString *s = new TString(Form("%d", i - 1));
    h1->Fill(s->Data(), realtime[i - 1]);
    h1->GetXaxis()->SetBinLabel(i, s->Data());
  }
  h1->Draw("b");
}
//______________________________________________________________________________
concurrent_queue::concurrent_queue(bool counter)
    : the_queue(), the_mutex(), the_condition_variable(&the_mutex), the_counter(0), nobjects(0),
      npriority(0) {
  // Concurrent queue constructor.
  the_counter = new TimeCounter(counter);
}

//______________________________________________________________________________
concurrent_queue::~concurrent_queue() {
  // destructor
  delete the_counter;
}

//______________________________________________________________________________
void concurrent_queue::push(TObject *data, bool priority) {
  // Push front and pop back policy for normal baskets, push back for priority
  // baskets.
  the_mutex.Lock();
  nobjects++;
  if (priority) {
    the_queue.push_back(data);
    npriority++;
  } else
    the_queue.push_front(data);
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
int concurrent_queue::size() const {
  the_mutex.Lock();
  int the_size = the_queue.size();
  the_mutex.UnLock();
  return the_size;
}

//______________________________________________________________________________
TObject *concurrent_queue::wait_and_pop() {
  //   --(*the_counter);
  the_mutex.Lock();
  while (the_queue.empty()) {
    //         Printf("WAITING");
    the_condition_variable.Wait();
  }

  TObject *popped_value = the_queue.back();
  the_queue.pop_back();
  nobjects--;
  if (npriority)
    npriority--;
  //   Printf("Popped basket %s", popped_value->GetName());
  //   ++(*the_counter);
  the_mutex.UnLock();
  return popped_value;
}

//______________________________________________________________________________
TObject *concurrent_queue::wait_and_pop_max(unsigned int nmax, unsigned int &n, TObject **array) {
  // Pop many objects in one go, maximum nmax. Will return the number of requested
  // objects, or the number available. The wait time for the client is minimal (at
  // least one object in the queue).
  the_mutex.Lock();
  while (the_queue.empty()) {
    the_condition_variable.Wait();
  }
  n = the_queue.size();
  if (n > nmax)
    n = nmax;
  unsigned int npopped = 0;
  while (npopped < n) {
    array[npopped++] = the_queue.back();
    the_queue.pop_back();
    nobjects--;
    if (npriority)
      npriority--;
  }
  the_mutex.UnLock();
  return array[0];
}

//______________________________________________________________________________
void concurrent_queue::pop_many(unsigned int n, TObject **array) {
  // Pop many objects in one go. Will keep the asking thread waiting until there
  // are enough objects in the queue.
  the_mutex.Lock();
  while (the_queue.size() < n) {
    the_condition_variable.Wait();
  }
  unsigned int npopped = 0;
  while (npopped < n) {
    array[npopped++] = the_queue.back();
    the_queue.pop_back();
    nobjects--;
    if (npriority)
      npriority--;
  }
  the_mutex.UnLock();
}

//______________________________________________________________________________
void concurrent_queue::Print() {
  if (the_counter)
    the_counter->Print();
}
#endif
