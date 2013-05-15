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
    Bool_t              AreOverlapping_v(const TGeoBBox_v  **box1, const TGeoMatrix  **mat1, const TGeoBBox_v  **box2, const TGeoMatrix  **mat2, Bool_t *  isin, const Int_t np);
    
   virtual void          ComputeNormal(Double_t *point, Double_t *dir, Double_t *norm);
   virtual void          ComputeNormal_v(Double_t *point, Double_t  *dir, Double_t  *norm, const Int_t np);
   virtual Bool_t        Contains(Double_t *point) const;
   virtual void          Contains_v(const Double_t *point, Bool_t *isin, Int_t np) const;
   static  Bool_t        Contains(const Double_t *point, Double_t dx, Double_t dy, Double_t dz, const Double_t *origin);
   virtual Bool_t        CouldBeCrossed(Double_t *point, Double_t *dir) const;
   virtual void          CouldBeCrossed_v(Double_t  *point, Double_t  *dir,Bool_t *  isin, const Int_t np) const;
    
   virtual Int_t         DistancetoPrimitive(Int_t px, Int_t py);

   virtual Double_t      DistFromInside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;

   inline  static  Double_t      DistFromInside(const Double_t *point,const Double_t *dir, 
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());

   static  void      DistFromInside_l(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax, Double_t * dist, int npoints);
   static  void      DistFromInside_v(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax, Double_t * dist, int npoints);


   virtual Double_t      DistFromOutside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;
   
   virtual void DistFromOutside_l(Double_t  *point, Double_t   *dir, Double_t *  isin, const Int_t np ) const;
   
   virtual void          DistFromOutside_v(Double_t   *point, Double_t  *dir, Int_t   *iact, Double_t    *step, Double_t    *safe, Double_t *  isin, const Int_t np ) const;
   static  Double_t      DistFromOutside(const Double_t *point,const Double_t *dir,
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());
   virtual Bool_t        GetPointsOnFacet(Int_t index, Int_t npoints, Double_t *array) const;
   virtual Bool_t        GetPointsOnSegments(Int_t npoints, Double_t *array) const;

   virtual Double_t      Safety(const Double_t *point, Bool_t in=kTRUE) const;
   void                  Safety_v(const Double_t *point, Double_t *safety, Int_t np, Bool_t in=kTRUE) const; 

   void                  Safety_v(const Double_t *point, Double_t *safety, const Bool_t * in, const Int_t np ) const; 

   ClassDef(TGeoBBox_v, 1)         // box primitive
};


inline
Double_t TGeoBBox_v::DistFromInside(const Double_t *point,const Double_t *dir,  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t /*stepmax*/)
{
// Computes distance from inside point to surface of the box.
// Boundary safe algorithm.

   Double_t s,smin,saf[6];
   Double_t newpt[3];
   Int_t i;

   newpt[0] = point[0] - origin[0];
   saf[0] = dx+newpt[0];
   saf[1] = dx-newpt[0];
   newpt[1] = point[1] - origin[1];
   saf[2] = dy+newpt[1];
   saf[3] = dy-newpt[1];
   newpt[2] = point[2] - origin[2];
   saf[4] = dz+newpt[2];
   saf[5] = dz-newpt[2];
   smin=TGeoShape::Big();

   double s1, s2;
   if (dir[0]!=0) 
     {
       // sign will be zero for negative numbers
       s = (dir[0]>0)? (saf[1]/dir[0]):(-saf[0]/dir[0]);
       if (s < smin) smin = s;
       if (s < 0) smin = 0.0;
     }


   if (dir[1]!=0) 
     {
       s = (dir[1]>0)? (saf[3]/dir[1]):(-saf[2]/dir[1]);
       if (s < smin) smin = s;
       if (s < 0) smin= 0.0;
     }   
   
   if (dir[2]!=0) 
     {
       s = (dir[2]>0)?(saf[5]/dir[2]):(-saf[4]/dir[2]);
       if (s < smin) smin = s;
       if (s < 0) smin = 0.0;
     }
    
    return smin;
}

#endif
