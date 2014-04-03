#ifndef GXTRACKLIASON_H
#define GXTRACKLIASON_H 

#include "GPTypeDef.h"
#include "GPMaterial.h"

struct
#ifdef __CUDA_
__align__(16)
#endif
GXTrackLiason
{
  GPMaterial *material; //material pts is persistent
};

#endif
