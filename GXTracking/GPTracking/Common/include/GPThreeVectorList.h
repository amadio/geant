#ifndef GPThreeVectorList_HH
#define GPThreeVectorList_HH 1

#include "GPThreeVector.h"

struct GPThreeVectorList
{
  int fSize;
  GPThreeVector vlist[8];
};

extern "C" {

FQUALIFIER 
void GPThreeVectorList_Constructor(GPThreeVectorList *This);

FQUALIFIER 
double GPThreeVectorList_size(GPThreeVectorList *This);

FQUALIFIER 
void GPThreeVectorList_pushback(GPThreeVectorList *This,
                                GPThreeVector v);

FQUALIFIER 
GPThreeVector GPThreeVectorList_popback(GPThreeVectorList *This, 
					int idx);

}

#endif
