// @(#)root/geom:$Id: TGeoBBox_v.h 27731 2009-03-09 17:40:56Z brun $
// Author: Andrei Gheata   24/10/01
   
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoBBox_v
#define ROOT_TGeoBBox_v

#ifndef ROOT_TGeoShape
#include "TGeoShape.h"
#endif

#include "TGeoBBox.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TGeoBBox_v - box class. All shape primitives inherit from this, their     //
//   constructor filling automatically the parameters of the box that bounds //
//   the given shape. Defined by 6 parameters :                              //
//      fDX, fDY, fDZ - half lengths on X, Y and Z axis                      //
//      fOrigin[3]    - position of box origin                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//#include <malloc.h>
#ifdef __INTEL_COMPILER
#include <immintrin.h> 
#endif

// this is used to pass data as struct of array
struct StructOfCoord
{
public:
  int np; // the size
  double *x;
  double *y;
  double *z;

  void fill(double * onedvec)
  {
    for(unsigned int i=0;i<np;++i)
      {
	x[i]=onedvec[3*i];
	y[i]=onedvec[3*i+1];
	z[i]=onedvec[3*i+2];
      }
  }

  void alloc(int np)
  {
    this->np=np;
#ifdef __INTEL_COMPILER
    x=(double*)_mm_malloc(sizeof(double)*np,32);
    y=(double*)_mm_malloc(sizeof(double)*np,32);
    z=(double*)_mm_malloc(sizeof(double)*np,32);
#else
    x=new double[np];
    y=new double[np];
    z=new double[np];
#endif
  }

  void dealloc()
  {
#ifdef __INTEL_COMPILER
    _mm_free(x);
    _mm_free(y);
    _mm_free(z);
#else
    delete[] x;
    delete[] y;
    delete[] z;
#endif
  }
};

class TGeoBBox_v : public TGeoBBox
{
public:
   // constructors
   TGeoBBox_v();
   TGeoBBox_v(Double_t dx, Double_t dy, Double_t dz, Double_t *origin=0);
   TGeoBBox_v(const char *name, Double_t dx, Double_t dy, Double_t dz, Double_t *origin=0);
   TGeoBBox_v(Double_t *param);
   // destructor
   virtual ~TGeoBBox_v();
   // methods
   static  Bool_t        AreOverlapping(const TGeoBBox_v *box1, const TGeoMatrix *mat1, const TGeoBBox_v *box2, const TGeoMatrix *mat2);
    
   virtual void          ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const;
   virtual void          ComputeNormal_v(const Double_t *point, const Double_t  *dir, Double_t  *norm, Int_t np) const;
   virtual void          ComputeNormal_l(const Double_t *point, const Double_t  *dir, Double_t  *norm, Int_t np) const;

   virtual Bool_t        Contains(const Double_t *point) const;
   virtual void          Contains_v(const Double_t *point, Bool_t *isin, Int_t np) const;
   virtual void          Contains_v(const StructOfCoord & point, Bool_t *isin, Int_t np) const;
   virtual void          Contains_l(const Double_t *point, Bool_t *isin, Int_t np) const;

   static  Bool_t        Contains(const Double_t *point, Double_t dx, Double_t dy, Double_t dz, const Double_t *origin);

   virtual Bool_t        CouldBeCrossed(const Double_t *point, const Double_t *dir) const;
   virtual void          CouldBeCrossed_l(const Double_t *point, const Double_t *dir, Bool_t * crossed, Int_t np ) const;
   virtual void          CouldBeCrossed_v(const Double_t  *point, const Double_t  *dir,Bool_t *  crossed, Int_t np) const;
   // SOA version
   virtual void          CouldBeCrossed_v(const StructOfCoord  & point, const StructOfCoord  & dir,Bool_t *  crossed, Int_t np) const;
    
   virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);

   virtual Double_t      DistFromInside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;

   static  Double_t      DistFromInside(const Double_t *point,const Double_t *dir, 
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());

   static  void      DistFromInside_l(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax, Double_t * dist, int npoints);
   static  void      DistFromInside_v(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax, Double_t * dist, int npoints);
   //SOA version
   static  void      DistFromInside_v(const StructOfCoord & point,const StructOfCoord & dir, 
   					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax, Double_t * dist, int npoints);



   virtual Double_t    DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   
   virtual void DistFromOutside_l(const Double_t  *point, const Double_t  *dir, Double_t *  distance, Int_t np ) const;
   
   virtual void          DistFromOutside_v(const Double_t *point, const Double_t  *dir, Int_t *iact, const Double_t  *step, Double_t  *safe, Double_t *  distance, Int_t np ) const;
 
 
   static  Double_t  DistFromOutside(const Double_t *point,const Double_t *dir,
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());
    // SOA version for static method
   static  void      DistFromOutside_v(const StructOfCoord & point,const StructOfCoord & dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
   // SOA trivial loop version for static method
   static  void      DistFromOutside_l(const StructOfCoord & point,const StructOfCoord & dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
   // trivial loop version for static method, no SOA
   static  void      DistFromOutside_l(const double * point,const double * dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);


   virtual Bool_t        GetPointsOnFacet(Int_t index, Int_t npoints, Double_t *array) const;
   virtual Bool_t        GetPointsOnSegments(Int_t npoints, Double_t *array) const;


   virtual Double_t      Safety(const Double_t *point, Bool_t in=kTRUE) const;
   virtual void          Safety_l(const Double_t *point, Double_t *safety, Int_t np, Bool_t in=kTRUE) const;

   virtual void          Safety_v(const Double_t *point, Double_t *safety, Int_t np, Bool_t in=kTRUE) const; 
   //SOA version
   virtual void          Safety_v(const StructOfCoord &point, Double_t *safety, Int_t np, Bool_t in=kTRUE) const;         

   virtual void          Safety_v(const Double_t *point, Double_t *safety, const Bool_t * in, const Int_t np ) const; 

   ClassDef(TGeoBBox_v, 1)         // box primitive
};



#endif
