#ifndef TEST_MIC_H
#define TEST_MIC_H

#include "GPTypeDef.h"
#include "GXFieldMap.h"

extern GLOBALFUNC void rk4_mic(GXFieldMap *magMap,
			       GXTrack *track,
			       unsigned int ntracks); 

extern GLOBALFUNC void rkf45_mic(GXFieldMap *magMap,
				 GXTrack *track,
				 unsigned int ntracks); 

extern GLOBALFUNC void nrk4_mic(GXFieldMap *magMap,
				GXTrack *track,
				unsigned int ntracks); 

extern GLOBALFUNC void rk4_cpu(GXFieldMap *magMap,
			       GXTrack *track,
			       unsigned int ntracks); 

extern GLOBALFUNC void rkf45_cpu(GXFieldMap *magMap,
				 GXTrack *track,
				 unsigned int ntracks); 

extern GLOBALFUNC void nrk4_cpu(GXFieldMap *magMap,
				GXTrack *track,
				unsigned int ntracks); 

#endif
