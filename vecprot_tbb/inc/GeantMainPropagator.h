#ifndef GEANT_MAINPROPAGATOR
#define GEANT_MAINPROPAGATOR

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantPropagator;

class GeantMainPropagator : public TObject
{
private:
	static GeantMainPropagator* fgInstance;
	static GeantPropagator* fgPropInstance;

public:
	GeantMainPropagator ();
	virtual ~GeantMainPropagator ();

	static GeantMainPropagator* Instance ();
	static GeantPropagator* PropInstance ();

	void SetParams (Int_t nthr, Int_t evtot, Int_t evbuf, Double_t tracksaver, Int_t maxperbask,
                  Int_t minFeeder, Int_t numPrior, Int_t dispThr,
                  Bool_t dbg=kFALSE, Int_t dbgTrk=-1);
	void Start (const char *geomfile="../geometry/cms.root", Bool_t graphics=kTRUE, Bool_t single=kFALSE);

	ClassDef(GeantMainPropagator,1)
};

#endif // GEANT_MAINPROPAGATOR

