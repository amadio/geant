
#include "UniformMagField.h"
// #include "PhysicalConstants.h"  //   For pi & twopi - Temporary solution ..
#include "SystemOfUnits.h"  //   For pi & twopi - Temporary solution ..

using geant::kPi;
using geant::kTwoPi;

// Virtual methods

void UniformMagField::ObtainFieldValue( const Vector3D<double> &position,
                                       Vector3D<double>  &fieldValue )
{
   GetFieldValue( position, fieldValue );
}

  /** @brief Vector interface for field retrieval */
void UniformMagField::ObtainFieldValueSIMD( const Vector3D<Double_v> &multiPosition, 
                                                  Vector3D<Double_v> &multiFieldValue )
{
   GetFieldValueSIMD( multiPosition, multiFieldValue );
}

// Constructor

UniformMagField::UniformMagField(double vField,
                                 double vTheta,
                                 double vPhi)
    :  VVectorField( 3, true )
{
   using namespace vecCore::math;
   if ( (vField<0) || (vTheta<0) || (vTheta>kPi) || (vPhi<0) || (vPhi>kTwoPi) )
   {
      // Exception("UniformMagField::UniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in UniformMagField::UniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << kPi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << kTwoPi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << kTwoPi << std::endl;
   }
   fFieldComponents.Set( vField*Sin(vTheta)*Cos(vPhi),
                         vField*Sin(vTheta)*Sin(vPhi),
                         vField*Cos(vTheta) );
}
