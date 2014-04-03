#ifndef UTIL_KERNEL_H
#define UTIL_KERNEL_H

#include "GXBrem.h"
#include "GXIoni.h"
#include "GXMsc.h"

extern "C" {

FQUALIFIER
void EMPhysics_Init(unsigned int tid,
		    curandState* devStates,
		    GXBrem *brem,
		    GXIoni *ioni,
		    GXMsc  *msc,
		    GPPhysicsTable* eBrem_table, 
		    GPPhysicsTable* eIoni_table, 
		    GPPhysicsTable* msc_table); 

FQUALIFIER
G4double EMPhysics_DefineStepLength(GXBrem *brem,
				    GXIoni *ioni,
				    GXMsc  *msc,
				    GXTrack *atrack,
				    GPForceCondition *condition,
				    GPGPILSelection  *selection); 
}

#endif
