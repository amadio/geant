#include "GeantVApplication.h"

//______________________________________________________________________________
GeantVApplication::GeantVApplication(GeantPropagator *prop):fPropagator(prop)
{
  // Ctor..

}

void GeantVApplication::setPropagator(GeantPropagator *prop){
	fPropagator=prop;
}
