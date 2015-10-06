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

	void SetParams (int nthr, int evtot, int evbuf, double tracksaver, int maxperbask,
                  int minFeeder, int numPrior, int dispThr,
                  bool dbg=false, int dbgTrk=-1);
	void Start (const char *geomfile="../geometry/cms.root", bool graphics=true, bool single=false);

	ClassDef(GeantMainPropagator,1)
};

#endif // GEANT_MAINPROPAGATOR

