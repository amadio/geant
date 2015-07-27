#include "GeantMainPropagator.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"


using std::max;

ClassImp (GeantMainPropagator)

GeantMainPropagator* GeantMainPropagator::fgInstance = 0;
GeantPropagator* GeantMainPropagator::fgPropInstance = 0;

GeantMainPropagator::GeantMainPropagator () :
	TObject ()
{
}

GeantMainPropagator::~GeantMainPropagator ()
{
}

GeantMainPropagator* GeantMainPropagator::Instance ()
{
	if (!fgInstance) fgInstance = new GeantMainPropagator();
	return fgInstance;
}

GeantPropagator* GeantMainPropagator::PropInstance ()
{
	if (!fgPropInstance) fgPropInstance = GeantPropagator::Instance();
	return fgPropInstance;
}

void GeantMainPropagator::SetParams (Int_t nthr, Int_t evtot, Int_t evbuf, Double_t tracksaver, Int_t maxperbask,
                                    Int_t minFeeder, Int_t numPrior, Int_t dispThr, Bool_t dbg, Int_t dbgTrk)
{
	GeantPropagator* p = GeantMainPropagator::PropInstance ();
   p->fNthreads = nthr;
	p->fNtotal = evtot;
	p->fNevents = evbuf;
	p->fNaverage = tracksaver;
	p->fNperBasket = maxperbask;

   p->fMinFeeder = max<int>(minFeeder, 2*nthr);
   p->fNevToPrioritize = numPrior;
   p->fDispThr = dispThr;

   p->fUseDebug = dbg;
   p->fDebugTrk = dbgTrk;

   WorkloadManager::Instance();
}

void GeantMainPropagator::Start (const char *geomfile, Bool_t graphics, Bool_t single)
{
   GeantPropagator* p = GeantMainPropagator::PropInstance ();
   p->PropagatorGeom(geomfile, graphics, single);
}
