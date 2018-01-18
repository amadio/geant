
#include "UniformMagField.h"
#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

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
   if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
   {
      // Exception("UniformMagField::UniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in UniformMagField::UniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << Constants::twopi << std::endl;
   }
   fFieldComponents.Set( vField*Sin(vTheta)*Cos(vPhi),
                         vField*Sin(vTheta)*Sin(vPhi),
                         vField*Cos(vTheta) );
}
