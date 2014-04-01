#ifndef GPGPU_KERNEL_H
#define GPGPU_KERNEL_H

#include "GPVPhysicalVolume.h"
#include "GPFieldMap.h"
#include "GPTrack.h"
#include "GPPhysicsTable.h"

// tracking_kernel wrapper
extern void tracking_gpu(GPVPhysicalVolume *world,
		      GPFieldMap *magMap,
		      GPPhysicsTable* eBrem_table,
		      GPPhysicsTable* eIoni_table,
		      GPPhysicsTable* msc_table,
                      GPTrack *track, 
		      size_t nTrackSize,
		      int nBlocks,
		      int nThreads); 

extern void tracking_cpu(GPVPhysicalVolume *world,
		      GPFieldMap *magMap,
		      GPPhysicsTable* eBrem_table,
		      GPPhysicsTable* eIoni_table,
		      GPPhysicsTable* msc_table,
		      GPTrack *trackIn,
		      size_t nTrackSize);

#endif
