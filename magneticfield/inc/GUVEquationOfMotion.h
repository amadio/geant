//
// class GUVEquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis     Dec 2014/Jan 2015
// -------------------------------------------------------------------

#ifndef GUV_EquationOfMotion_H
#define GUV_EquationOfMotion_H

#include <cassert>
#include <iostream>

// #include <vector>
#include "base/Vector3D.h"
#include <Geant/Config.h>  // To define GEANT_FORCE_INLINE

// #include "GUVTypes.hh"      // "globals.hh"
#include "GUVField.h"   // required in inline method implementations

class GUVEquationOfMotion 
{
  public:  // with description

     GUVEquationOfMotion( GUVField *Field, unsigned short verbose=0 );
     virtual ~GUVEquationOfMotion();
       // Constructor and virtual destructor. No operations, just checks

     virtual void EvaluateRhsGivenB( const  double     yVec[],
                                     const  vecgeom::Vector3D<double> &B,
                                            double     charge,
                                            double     dydx[] ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

     // virtual void SetChargeMomentumMass(double particleCharge,
     //                                    double MomentumXc,
     //                                    double MassXc2) = 0;
     //   // Set the charge, momentum and mass of the current particle
     //   // --> used to set the equation's coefficients ...

     inline void RightHandSide( const  double y[],
                                        double charge ,
                                        double dydx[] ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // It uses the virtual function EvaluateRhsGivenB

     inline void EvaluateRhsReturnB( const double y[],
                              double       dydx[],
                              double       charge,
                  vecgeom::Vector3D<double> &field ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     inline
     void GetFieldValue( const vecgeom::Vector3D<double> &position,
                         vecgeom::Vector3D<double>       &fieldValue ) const;

     const GUVField* GetFieldObj() const {return fField;}
           GUVField* GetFieldObj()       {return fField;}
     void            SetFieldObj(GUVField* pField){fField=pField;}

     bool         Initialised() const { return fInitialised; } 
     unsigned int GetId() const       { return fEquationId; }
     static unsigned int GetNumCreated() { return fNumObjectsCreated; }
     static unsigned int GetNumLive() { return fNumObjectsCreated - fNumObjectsDeleted; }
       // For debugging, checking

     friend std::ostream&
             operator<<( std::ostream& os, const GUVEquationOfMotion& eq);

  public:
     static const unsigned int idxTime=3;  // Convention for location of time 't' in vector

  private:
     static unsigned int fNumObjectsCreated;
     static unsigned int fNumObjectsDeleted;
     // const int GUVmaximum_number_of_field_components = 24;
     enum { GUVmaximum_number_of_field_components = 24 } ;

     GUVField *     fField;
     unsigned int   fEquationId;  //
     unsigned short fVerbose;
     bool           fInitialised;
};

// #include "GUVEquationOfMotion.icc"

//  Inline implementation
//
// -------------------------------------------------------------------

inline
GUVEquationOfMotion::GUVEquationOfMotion(GUVField* pField, unsigned short verbose)
   : fField(pField), fEquationId(fNumObjectsCreated++),
     fVerbose(verbose), fInitialised(false)
{
   if( fVerbose )
      std::cout << " Created Equation " << this << " info= " << *this << std::endl;
}

inline
void GUVEquationOfMotion::GetFieldValue( const vecgeom::Vector3D<double> &Position,
                                               vecgeom::Vector3D<double>  &FieldValue ) const
{
   fField-> GetFieldValue( Position, FieldValue );
}

GEANT_FORCE_INLINE
void
GUVEquationOfMotion::RightHandSide( const  double y[],
                                           double charge,
                                           double dydx[]  ) const
{
   using ThreeVectorD = vecgeom::Vector3D<double>;

   // double Field[GUVmaximum_number_of_field_components];
   ThreeVectorD  field;
   // double PositionAndTime[4];

   ThreeVectorD  position( y[0], y[1], y[2] );

   //  PositionAndTime[0] = y[0];
   //  PositionAndTime[1] = y[1];
   //  PositionAndTime[2] = y[2];
   // Global Time -- ignored for now
   //  PositionAndTime[3] = y[idxTime];  // See GUVFieldTrack::LoadFromArray

   GetFieldValue( position, field );
   // GetFieldValue( y, Field_3vf );   
   EvaluateRhsGivenB( y, field, charge, dydx );
}

GEANT_FORCE_INLINE
void
GUVEquationOfMotion::
EvaluateRhsReturnB( const double           y[],
                          double          dydx[],
                          double          charge,
                    vecgeom::Vector3D<double> &field
                  ) const
{
   using ThreeVector = vecgeom::Vector3D<double>;
   
   GetFieldValue( ThreeVector(y[0], y[1], y[2]), field) ;
   EvaluateRhsGivenB( y, field, charge, dydx );
}

#endif /* GUV_EquationOfMotion_DEF */
