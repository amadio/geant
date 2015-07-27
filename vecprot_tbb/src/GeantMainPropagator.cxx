#include "GeantMainPropagator.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"

using std::max;

ClassImp(GeantMainPropagator)

    GeantMainPropagator *GeantMainPropagator::fgInstance = 0;
GeantPropagator *GeantMainPropagator::fgPropInstance = 0;

GeantMainPropagator::GeantMainPropagator() : TObject() {}

GeantMainPropagator::~GeantMainPropagator() {}

GeantMainPropagator *GeantMainPropagator::Instance() {
  if (!fgInstance)
    fgInstance = new GeantMainPropagator();
  return fgInstance;
}

GeantPropagator *GeantMainPropagator::PropInstance() {
  if (!fgPropInstance)
    fgPropInstance = GeantPropagator::Instance();
  return fgPropInstance;
}

void GeantMainPropagator::SetParams(int nthr, int evtot, int evbuf, double tracksaver, int maxperbask,
                                    int minFeeder, int numPrior, int dispThr, Bool_t dbg, int dbgTrk) {
  GeantPropagator *p = GeantMainPropagator::PropInstance();
  p->fNthreads = nthr;
  p->fNtotal = evtot;
  p->fNevents = evbuf;
  p->fNaverage = tracksaver;
  p->fNperBasket = maxperbask;

  p->fMinFeeder = max<int>(minFeeder, 2 * nthr);
  p->fNevToPrioritize = numPrior;
  p->fDispThr = dispThr;

  p->fUseDebug = dbg;
  p->fDebugTrk = dbgTrk;

  WorkloadManager::Instance();
}

void GeantMainPropagator::Start(const char *geomfile, Bool_t graphics, Bool_t single) {
  GeantPropagator *p = GeantMainPropagator::PropInstance();
  p->PropagatorGeom(geomfile, graphics, single);
}
