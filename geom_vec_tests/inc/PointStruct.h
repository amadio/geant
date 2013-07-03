#ifndef POINTSTRUCT_H
#define POINTSTRUCT_H

#ifdef __INTEL_COMPILER
#include <immintrin.h> 
#else
#include "mm_malloc.h"
#endif


///////////////////////////////////////////////////////////////////
//
//  Class to store vectors of 3D dimensional (vectors) data in SOA
//  form
//
///////////////////////////////////////////////////////////////////


// this is used to pass data as struct of array
struct StructOfCoord
{
private:
  double * xvec;
  double * yvec;
  double * zvec;

public:
  int np; // the size
  double *x; // these are just "views" to the real data in xvec ( with the idea that we can set x to point to different start locations in xvec )
  double *y;
  double *z;

  void setstartindex(int index)
  {
    x=&xvec[index];
    y=&yvec[index];
    z=&zvec[index];
  }
  
  // checks if index will be aligned and if not go back to an aligned address
  void setstartindex_aligned(int index)
  {
    // TOBEDONE
  }

  void fill(double * onedvec)
  {
    for(unsigned int i=0;i<np;++i)
      {
	xvec[i]=onedvec[3*i];
	yvec[i]=onedvec[3*i+1];
	zvec[i]=onedvec[3*i+2];
      }
  }

  void alloc(int np)
  {
    this->np=np;
    xvec=(double*)_mm_malloc(sizeof(double)*np,32); // aligned malloc (32 for AVX )
    yvec=(double*)_mm_malloc(sizeof(double)*np,32);
    zvec=(double*)_mm_malloc(sizeof(double)*np,32);
    x=xvec;y=yvec;z=zvec;
  }

  void dealloc()
  {
    _mm_free(xvec);
    _mm_free(yvec);
    _mm_free(zvec);
    x=0;y=0;z=0;
  }
};

#endif // POINTSTRUCT_H
