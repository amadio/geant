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
   TGeoBBox_v(double dx, double dy, double dz, double *origin=0);
   TGeoBBox_v(const char *name, double dx, double dy, double dz, double *origin=0);
   TGeoBBox_v(double *param);
   virtual ~TGeoBBox_v();    // destructor

   // methods
   static  bool        AreOverlapping(const TGeoBBox_v *box1, const TGeoMatrix *mat1, const TGeoBBox_v *box2, const TGeoMatrix *mat2);
    
   virtual void          ComputeNormal(const double *point, const double *dir, double *norm) const;
   virtual void          ComputeNormal_v(const double *point, const double  *dir, double  *norm, int np) const;
   virtual void          ComputeNormal_l(const double *point, const double  *dir, double  *norm, int np) const;

   virtual bool        Contains(const double *point) const;
   virtual void          Contains_v(const double *point, bool *isin, int np) const;
   virtual void          Contains_v(const StructOfCoord & point, bool *isin, int np) const;
   virtual void          Contains_l(const double *point, bool *isin, int np) const;

   static  bool        Contains(const double *point, double dx, double dy, double dz, const double *origin);

   virtual bool        CouldBeCrossed(const double *point, const double *dir) const;
   virtual void          CouldBeCrossed_l(const double *point, const double *dir, bool * crossed, int np ) const;
   virtual void          CouldBeCrossed_v(const double  *point, const double  *dir,bool *  crossed, int np) const;
   // SOA version
   virtual void          CouldBeCrossed_v(const StructOfCoord  & point, const StructOfCoord  & dir,bool *  crossed, int np) const;
    
   virtual int         DistancetoPrimitive(int px, int py);

   // DISTFROMINSIDE
   virtual double      DistFromInside(double *point, double *dir, int iact=1, 
                                   double step=TGeoShape::Big(), double *safe=0) const;

   // often this method will be the general entry point from outside but we just dispatch it to the static method and ignore iact and safe )
   virtual void    DistFromInside_v(StructOfCoord const & point, StructOfCoord const & dir, int iact, 
				     double const * step, double *safe, double * distance , int np) const {
     TGeoBBox_v::DistFromInsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(), step, distance, np );
   }


   static  double  DistFromInsideS(const double *point,const double *dir, 
                                   double dx, double dy, double dz, const double *origin, double stepmax=TGeoShape::Big());

   static  void      DistFromInsideS_l(const double *point,const double *dir, 
					  double dx, double dy, double dz, const double *origin, const double *stepmax, double * dist, int npoints);
   // SOA version for _l
   static  void      DistFromInsideS_l(const StructOfCoord & point,const StructOfCoord & dir, 
					  double dx, double dy, double dz, const double *origin, const double *stepmax, double * dist, int npoints);


   static  void      DistFromInsideS_v(const double *point,const double *dir, 
					  double dx, double dy, double dz, const double *origin, const double *stepmax, double * dist, int npoints);
   //SOA version
   static  void      DistFromInsideS_v(const StructOfCoord & point,const StructOfCoord & dir, 
   					  double dx, double dy, double dz, const double *origin, const double *stepmax, double * dist, int npoints);


   // DISTFROMOUTSIDE
   virtual double    DistFromOutside(const double *point, const double *dir, int iact=1, 
                                   double step=TGeoShape::Big(), double *safe=0) const;

   // often this method will be the general entry point from outside ( we dispatch it to the static method )
   virtual void    DistFromOutside_v(StructOfCoord const & point, StructOfCoord const & dir, int iact, 
				     double const * step, double *safe, double * distance , int np) const {
     TGeoBBox_v::DistFromOutsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
   }

#ifdef VEC_EXTENSIONS
   virtual void    DistFromOutside_Masked_v(double const *m, StructOfCoord const & point, StructOfCoord const & dir, int iact, 
				     double const * step, double *safe, double * distance , int np) const {
     TGeoBBox_v::DistFromOutsideS_Masked_v(m, point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
   }
#endif

   // static version ( added "S" to it to make it consistent with other TGeo classes )
   static  double  DistFromOutsideS(const double *point,const double *dir,
                                   double dx, double dy, double dz, const double *origin, double stepmax=TGeoShape::Big());
   
   virtual void DistFromOutside_l(const double  *point, const double  *dir, double *  distance, int np ) const;
 
    // SOA version for static method
   static  void      DistFromOutsideS_v(const StructOfCoord & point,const StructOfCoord & dir,
					   double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);
  // SOA version for static method
#ifdef VEC_EXTENSIONS
   static  void      DistFromOutsideS_Masked_v(double const*, const StructOfCoord & point,const StructOfCoord & dir,
					   double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);
#endif   

   // SOA trivial loop version for static method
   static  void      DistFromOutsideS_l(const StructOfCoord & point,const StructOfCoord & dir,
					   double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);
   // trivial loop version for static method, no SOA
   static  void      DistFromOutsideS_l(const double * point,const double * dir,
					   double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);

   // this code is meant to try a vectorization approach at a lower level with vector extensions or intrinsics  
#ifdef VEC_EXTENSIONS
   // public:
   //   static void      DistFromOutsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
   //
   //                     double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);
 private:
   // this is the actual kernel doing the computation with possibility of early return
   static void DistFromOutsideS_v4( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance );
  
 public:
#endif




   virtual bool        GetPointsOnFacet(int index, int npoints, double *array) const;
   virtual bool        GetPointsOnSegments(int npoints, double *array) const;


   virtual double      Safety(const double *point, bool in=true) const;
   virtual void          Safety_l(const double *point, double *safety, int np, bool in=true) const;

   virtual void          Safety_v(const double *point, bool in, double *safety, int np ) const; 
   //SOA version
   virtual void          Safety_v(const StructOfCoord &point, bool in, double *safety, int np ) const;         

   virtual void          Safety_v(const double *point, double *safety, const bool * in, const int np ) const; 

   //   ClassDef(TGeoBBox_v, 1)         // box primitive
};



#endif
