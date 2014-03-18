//A minimal static implementation for GPBox and GPGenericTrap 
//that have 8 vertices

#include "GPThreeVectorList.h"

FQUALIFIER 
void GPThreeVectorList_Constructor(GPThreeVectorList *This)
{
  This->fSize = 0;
}

FQUALIFIER 
double GPThreeVectorList_size(GPThreeVectorList *This)
{
  return This->fSize;
}

FQUALIFIER 
void GPThreeVectorList_pushback(GPThreeVectorList *This,
				GPThreeVector v)
{
  if(This->fSize < 8) This->vlist[This->fSize] = v;  
  (This->fSize)++;
}

FQUALIFIER 
GPThreeVector GPThreeVectorList_popback(GPThreeVectorList *This, int idx)
{
  if(idx < 8) return This->vlist[idx];  
  else return GPThreeVector_create(0,0,0);
}
