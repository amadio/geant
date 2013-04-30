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
// TGeoBBox_v - box class. All shape primitives inherit from this, their       //
//   constructor filling automatically the parameters of the box that bounds //
//   the given shape. Defined by 6 parameters :                              //
//      fDX, fDY, fDZ - half lengths on X, Y and Z axis                      //
//      fOrigin[3]    - position of box origin                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

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
   virtual void          ComputeNormal(Double_t *point, Double_t *dir, Double_t *norm);
   virtual Bool_t        Contains(Double_t *point) const;
   virtual void          Contains_v(const Double_t *point, Bool_t *isin, Int_t np) const;
   static  Bool_t        Contains(const Double_t *point, Double_t dx, Double_t dy, Double_t dz, const Double_t *origin);
   virtual Bool_t        CouldBeCrossed(Double_t *point, Double_t *dir) const;
   virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);
   virtual Double_t      DistFromInside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   static  Double_t      DistFromInside(const Double_t *point,const Double_t *dir, 
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());
   virtual Double_t      DistFromOutside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   static  Double_t      DistFromOutside(const Double_t *point,const Double_t *dir, 
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());
   virtual Bool_t        GetPointsOnFacet(Int_t index, Int_t npoints, Double_t *array) const;
   virtual Bool_t        GetPointsOnSegments(Int_t npoints, Double_t *array) const;
   virtual Double_t      Safety(Double_t *point, Bool_t in=kTRUE) const;

   ClassDef(TGeoBBox_v, 1)         // box primitive
};

#endif
