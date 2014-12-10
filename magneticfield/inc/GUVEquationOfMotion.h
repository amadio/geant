//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GUVEquationOfMotion.hh 71664 2013-06-20 08:36:05Z gcosmo $
//
//
// class GUVEquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis
// -------------------------------------------------------------------

#ifndef GUV_EquationOfMotion_DEF
#define GUV_EquationOfMotion_DEF

// #include "GUVTypes.hh"      // "globals.hh"
#include "GUVField.h"   // required in inline method implementations

class GUVEquationOfMotion 
{
  public:  // with description

     GUVEquationOfMotion( GUVField *Field );
     virtual ~GUVEquationOfMotion();
       // Constructor and virtual destructor. No operations.

     virtual void EvaluateRhsGivenB( const  double y[],
                                     const  double B[3],
                                     double dydx[] ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

     virtual void SetChargeMomentumMass(double particleCharge,
                                        double MomentumXc,
                                        double MassXc2) = 0;
       // Set the charge, momentum and mass of the current particle
       // --> used to set the equation's coefficients ...

     inline
     void RightHandSide( const  double y[],
                                double dydx[] ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // (It is not virtual, but calls the virtual function above.)

     void EvaluateRhsReturnB( const  double y[],
                              double dydx[],
                              double Field[]  ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     void GetFieldValue( const  double Point[4],
                                double Field[] )  const;
       // Obtain only the field - the stepper assumes it is pure Magnetic.
       // Not protected, because GUVRKG3_Stepper uses it directly.

     const GUVField* GetFieldObj() const;
     void           SetFieldObj(GUVField* pField);

  private:
     // const int GUVmaximum_number_of_field_components = 24;
     enum { GUVmaximum_number_of_field_components = 24 } ;

     GUVField *itsField;

};

#include "GUVEquationOfMotion.icc"

#endif /* GUV_EquationOfMotion_DEF */
