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
   TGeoBBox_v(double dx, double dy, double dz, double *origin=0);
   TGeoBBox_v(const char *name, double dx, double dy, double dz, double *origin=0);
   TGeoBBox_v(double *param);
   // destructor
   virtual ~TGeoBBox_v();
   // methods
   static  Bool_t        AreOverlapping(const TGeoBBox_v *box1, const TGeoMatrix *mat1, const TGeoBBox_v *box2, const TGeoMatrix *mat2);
   virtual void          ComputeNormal(double *point, double *dir, double *norm);
   virtual Bool_t        Contains(double *point) const;
   virtual void          Contains_v(const double *point, Bool_t *isin, Int_t np) const;
   static  Bool_t        Contains(const double *point, double dx, double dy, double dz, const double *origin);
   virtual Bool_t        CouldBeCrossed(double *point, double *dir) const;
   virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);
   virtual double      DistFromInside(double *point, double *dir, Int_t iact=1, 
                                   double step=TGeoShape::Big(), double *safe=0) const;
   static  double      DistFromInside(const double *point,const double *dir, 
                                   double dx, double dy, double dz, const double *origin, double stepmax=TGeoShape::Big());
   virtual double      DistFromOutside(double *point, double *dir, Int_t iact=1, 
                                   double step=TGeoShape::Big(), double *safe=0) const;
   static  double      DistFromOutside(const double *point,const double *dir, 
                                   double dx, double dy, double dz, const double *origin, double stepmax=TGeoShape::Big());
   virtual Bool_t        GetPointsOnFacet(Int_t index, Int_t npoints, double *array) const;
   virtual Bool_t        GetPointsOnSegments(Int_t npoints, double *array) const;
   virtual double      Safety(double *point, Bool_t in=kTRUE) const;

   ClassDef(TGeoBBox_v, 1)         // box primitive
};

#endif
