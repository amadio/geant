//
// class GUVVectorEquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis     Dec 2014/Jan 2015
// -------------------------------------------------------------------

#ifndef GUV_VectorEquationOfMotion_H
#define GUV_VectorEquationOfMotion_H

#include <cassert>
#include <iostream>

// #include <vector>
#include "backend/Backend.h"   // From VecGeom
#include "base/Vector3D.h"

// #include "GUVTypes.hh"      // "globals.hh"
#include "GUVField.h"   // required in inline method implementations
#include "GUVVectorField.h"
#include "GUVEquationOfMotion.h"

#include "VcFloatBackend.h"

class GUVVectorEquationOfMotion //: public GUVEquationOfMotion 
{
  
  typedef typename vecgeom::kVc::precision_v      Double_v;
  typedef typename vecgeom::kVcFloat::precision_v Float_v;
//  Tried alternative: 
//     typedef typename vecgeom::kVc::precision_v      kVcFloat::precision_v

  public:  // with description

     GUVVectorEquationOfMotion( GUVVectorField *Field, unsigned int verbose=0 );
     virtual ~GUVVectorEquationOfMotion();     
       // Constructor and virtual destructor. No operations, just checks

     virtual
     void EvaluateRhsGivenB( const  Double_v yVec[],
                             const  vecgeom::Vector3D<Float_v> B,  // Was double B[3],
                             const  Double_v charge, 
                                    Double_v dydx[]             ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

      
      virtual void InitializeCharge(double particleCharge)=0;
       // Must be called to correctly initialise and provide charge
      virtual void InvalidateParameters()=0;
      
      inline void InformReady(); // All parameters have been set (charge+)
      inline void InformDone();  // Invalidate charge, other parameters
      inline void CheckInitialization() const; // Ensure initialization
      inline void CheckDone() const;

     // virtual void SetChargeMomentumMass(double particleCharge,
     //                                    double MomentumXc,
     //                                    double MassXc2) = 0;
     //   // Set the charge, momentum and mass of the current particle
     //   // --> used to set the equation's coefficients ...

      inline void RightHandSide( const  Double_v y[],
                                        Double_v charge ,
                                        Double_v dydx[] ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // It uses the virtual function EvaluateRhsGivenB

     void EvaluateRhsReturnB( const Double_v y[],
                                    Double_v dydx[],
                                    Double_v charge,
                                    vecgeom::Vector3D<Float_v> &Field ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     void GetFieldValue( const Double_v Point[4],
                               Double_v Field[] )  const;
       // Obtain only the field - the stepper assumes it is pure Magnetic.
       // Not protected, because GUVRKG3_Stepper uses it directly.
     inline
     void GetFieldValue( const  Double_v                    Point[4],
                                vecgeom::Vector3D<Float_v>  &FieldValue ) const;

     inline
     void GetFieldValue( const vecgeom::Vector3D<Double_v> &Position,
                               vecgeom::Vector3D<Float_v>  &FieldValue ) const;

     const GUVVectorField* GetFieldObj() const {return fField;}
           GUVVectorField* GetFieldObj()       {return fField;}
     void                  SetFieldObj(GUVVectorField* pField){fField=pField;}

     bool         Initialised() const { return fInitialised; } 
     unsigned int GetId() const       { return fEquationId; }
     static unsigned int GetNumCreated() { return fNumObjectsCreated; }
     static unsigned int GetNumLive() { return fNumObjectsCreated - fNumObjectsDeleted; }


     friend std::ostream&
             operator<<( std::ostream& os, const GUVVectorEquationOfMotion& eq);


  public:
     static const unsigned int idxTime=3;  // Convention for location of time 't' in vector

  private:
     static unsigned int fNumObjectsCreated;
     static unsigned int fNumObjectsDeleted;
     // const int GUVmaximum_number_of_field_components = 24;
     enum { GUVmaximum_number_of_field_components = 24 } ;

     GUVVectorField *     fField;
     unsigned int   fEquationId;  //
     unsigned short fVerbose;
     bool           fInitialised;
};


inline
GUVVectorEquationOfMotion::GUVVectorEquationOfMotion(GUVVectorField* pField, unsigned int verbose)
   : // GUVEquationOfMotion(pField, verbose),
     fField(pField), 
     fEquationId(fNumObjectsCreated++),
     fVerbose(verbose), 
     fInitialised(false)
{
   if(fVerbose)
   {
     std::cout << " Called GUVVectorEquationOfMotion::InformDone() " << std::endl;
     std::cout << *this << std::endl;
   }
}

inline
void GUVVectorEquationOfMotion::InformReady() // was Initialize()
{
   fInitialised= true;
}

inline
void GUVVectorEquationOfMotion::InformDone()  // was Clear() and before Finished();
{
   if(fVerbose)
   {
     std::cout << " Called GUVVectorEquationOfMotion::InformDone() " << std::endl;
     std::cout << *this << std::endl;
   }
   assert( fInitialised );
   fInitialised= false;
}

inline
void GUVVectorEquationOfMotion::GetFieldValue( const typename vecgeom::kVc::precision_v Point[4],
                                                     typename vecgeom::kVc::precision_v Field[] ) const
{
   vecgeom::Vector3D<typename vecgeom::kVc::precision_v> Position( Point[0], Point[1], Point[2] );
   // vecgeom::Vector3D<Double_v> Position( Point[0], Point[1], Point[2] );   
   vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v>  FieldVec;
   // vecgeom::Vector3D<Float_v>  FieldVec;   
   fField-> GetFieldValue( Position, FieldVec );
   Field[0] = (Double_v) FieldVec[0];
   Field[1] = (Double_v) FieldVec[1];
   Field[2] = (Double_v) FieldVec[2];
}

inline
void GUVVectorEquationOfMotion::GetFieldValue( const Double_v Point[4], // Was: const vecgeom::Vector3D<double> &Position,
                                               // Tried alternative: vecgeom::Vector3D<Float_v>  &FieldValue
                                               vecgeom::Vector3D<vecgeom::kVcFloat::precision_v>  &FieldValue                                               
   ) const   
{
   vecgeom::Vector3D<Double_v> Position( Point[0], Point[1], Point[2] );
   fField-> GetFieldValue( Position, FieldValue );
}

inline
void GUVVectorEquationOfMotion::GetFieldValue( const vecgeom::Vector3D<typename vecgeom::kVc::precision_v>      &Position,
                                                     vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v> &FieldValue                                                          // Tried alternative: vecgeom::Vector3D<Float_v>  &FieldValue
   ) const
{
   fField-> GetFieldValue( Position, FieldValue );
}

inline
void
GUVVectorEquationOfMotion::RightHandSide( const typename vecgeom::kVc::precision_v y[],
                                                typename vecgeom::kVc::precision_v charge,
                                                typename vecgeom::kVc::precision_v dydx[] ) const
{
   using ThreeVectorF = vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v>;
   // using ThreeVectorF = vecgeom::Vector3D<Float_v>;  // Tried alternative
   using ThreeVectorD = vecgeom::Vector3D<Double_v>; // Was: typename vecgeom::kVc::precision_v>;
   CheckInitialization();

   ThreeVectorF  Field_3vf;
   ThreeVectorD  Position( y[0], y[1], y[2] );

   GetFieldValue( Position, Field_3vf );
   EvaluateRhsGivenB( y, Field_3vf, charge, dydx );
}

#include <iostream>

void GUVVectorEquationOfMotion::CheckInitialization() const
{
#ifdef GUVERBOSE
   if( fVerbose && !fInitialised ){
      std::cerr << "GUVVectorEquationOfMotion is not Initialised" << std::endl;
   }
#endif
   assert( fInitialised );
}

void GUVVectorEquationOfMotion::CheckDone() const
{
#ifdef GUVERBOSE
   if( fVerbose && fInitialised ){
      std::cerr << "GUVVectorEquationOfMotion was NOT told it is Done!" << std::endl;
   }
#endif
   assert( !fInitialised );
}


#endif /* GUV_VectorEquationOfMotion_DEF */
