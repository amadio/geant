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

#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"
#endif

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TGeoBBox_v - box class. All shape primitives inherit from this, their     //
//   constructor filling automatically the parameters of the box that bounds //
//   the given shape. Defined by 6 parameters :                              //
//      fDX, fDY, fDZ - half lengths on X, Y and Z axis                      //
//      fOrigin[3]    - position of box origin                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "PointStruct.h" // for SOA data

class TGeoBBox_v : public TGeoBBox
{
public:
   // constructors
   TGeoBBox_v();
   TGeoBBox_v(Double_t dx, Double_t dy, Double_t dz, Double_t *origin=0);
   TGeoBBox_v(const char *name, Double_t dx, Double_t dy, Double_t dz, Double_t *origin=0);
   TGeoBBox_v(Double_t *param);
   virtual ~TGeoBBox_v();    // destructor

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

   // DISTFROMINSIDE
   virtual Double_t      DistFromInside(Double_t *point, Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;

   // often this method will be the general entry point from outside but we just dispatch it to the static method and ignore iact and safe )
   virtual void    DistFromInside_v(StructOfCoord const & point, StructOfCoord const & dir, Int_t iact, 
				     Double_t const * step, Double_t *safe, Double_t * distance , Int_t np) const {
     TGeoBBox_v::DistFromInsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(), step, distance, np );
   }


   static  Double_t  DistFromInsideS(const Double_t *point,const Double_t *dir, 
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());

   static  void      DistFromInsideS_l(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t *stepmax, Double_t * dist, int npoints);
   // SOA version for _l
   static  void      DistFromInsideS_l(const StructOfCoord & point,const StructOfCoord & dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t *stepmax, Double_t * dist, int npoints);


   static  void      DistFromInsideS_v(const Double_t *point,const Double_t *dir, 
					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t *stepmax, Double_t * dist, int npoints);
   //SOA version
   static  void      DistFromInsideS_v(const StructOfCoord & point,const StructOfCoord & dir, 
   					  Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t *stepmax, Double_t * dist, int npoints);


   // DISTFROMOUTSIDE
   virtual Double_t    DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact=1, 
                                   Double_t step=TGeoShape::Big(), Double_t *safe=0) const;

   // often this method will be the general entry point from outside ( we dispatch it to the static method )
   virtual void    DistFromOutside_v(StructOfCoord const & point, StructOfCoord const & dir, Int_t iact, 
				     Double_t const * step, Double_t *safe, Double_t * distance , Int_t np) const {
     TGeoBBox_v::DistFromOutsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
   }

#ifdef VEC_EXTENSIONS
   virtual void    DistFromOutside_Masked_v(double const *m, StructOfCoord const & point, StructOfCoord const & dir, Int_t iact, 
				     Double_t const * step, Double_t *safe, Double_t * distance , Int_t np) const {
     TGeoBBox_v::DistFromOutsideS_Masked_v(m, point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
   }
#endif

   // static version ( added "S" to it to make it consistent with other TGeo classes )
   static  Double_t  DistFromOutsideS(const Double_t *point,const Double_t *dir,
                                   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, Double_t stepmax=TGeoShape::Big());
   
   virtual void DistFromOutside_l(const Double_t  *point, const Double_t  *dir, Double_t *  distance, Int_t np ) const;
 
    // SOA version for static method
   static  void      DistFromOutsideS_v(const StructOfCoord & point,const StructOfCoord & dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
  // SOA version for static method
#ifdef VEC_EXTENSIONS
   static  void      DistFromOutsideS_Masked_v(double const*, const StructOfCoord & point,const StructOfCoord & dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
#endif   

   // SOA trivial loop version for static method
   static  void      DistFromOutsideS_l(const StructOfCoord & point,const StructOfCoord & dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
   // trivial loop version for static method, no SOA
   static  void      DistFromOutsideS_l(const double * point,const double * dir,
					   Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);

   // this code is meant to try a vectorization approach at a lower level with vector extensions or intrinsics  
#ifdef VEC_EXTENSIONS
   // public:
   //   static void      DistFromOutsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
   //
   //                     Double_t dx, Double_t dy, Double_t dz, const Double_t *origin, const Double_t * stepmax, Double_t * distance, Int_t np);
 private:
   // this is the actual kernel doing the computation with possibility of early return
   static void DistFromOutsideS_v4( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance );
  
 public:
#endif




   virtual Bool_t        GetPointsOnFacet(Int_t index, Int_t npoints, Double_t *array) const;
   virtual Bool_t        GetPointsOnSegments(Int_t npoints, Double_t *array) const;


   virtual Double_t      Safety(const Double_t *point, Bool_t in=kTRUE) const;
   virtual void          Safety_l(const Double_t *point, Double_t *safety, Int_t np, Bool_t in=kTRUE) const;

   virtual void          Safety_v(const Double_t *point, Bool_t in, Double_t *safety, Int_t np ) const; 
   //SOA version
   virtual void          Safety_v(const StructOfCoord &point, Bool_t in, Double_t *safety, Int_t np ) const;         

   virtual void          Safety_v(const Double_t *point, Double_t *safety, const Bool_t * in, const Int_t np ) const; 

   //   ClassDef(TGeoBBox_v, 1)         // box primitive
};



#endif
