//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef UniformMagField_H
#define UniformMagField_H

#include <iostream>

#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

class UniformMagField
{
public:

  using Double_v = Geant::Double_v;
  using Float_v =  Geant::Float_v;
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  static constexpr int   gNumFieldComponents= 3;
  static constexpr bool  gFieldChangesEnergy= false;
  
  /** @brief Constructor providing the constant field value (cartesian) */
  UniformMagField( const vecgeom::Vector3D<float>& fieldVector )
    : fFieldComponents(fieldVector) {}

  /** @brief Constructor providing the constant field value (spherical) */
  UniformMagField( double vField,
                   double vTheta,
                   double vPhi  );

  /** @brief Destructor */
  ~UniformMagField() {}

  /** @brief Scalar interface for field retrieval */
  void  GetFieldValue( const Vector3D<double> &position, 
                             Vector3D<double> &fieldValue )
  {
    GetFieldValue<double>(position, fieldValue);
  }

  /** @brief Vector interface for field retrieval */
  void GetFieldValueSIMD( const Vector3D<Double_v> &position, 
                                Vector3D<Double_v> &fieldValue )
  {
    GetFieldValue<Double_v>(position, fieldValue);
  }

  /** @brief Templated field interface */
  template <typename Real_v>
  void GetFieldValue( const Vector3D<Real_v> & /*position*/,
                            Vector3D<Real_v> &fieldValue )
  {
    fieldValue.Set(Real_v(fFieldComponents.x()),
                   Real_v(fFieldComponents.y()),
                   Real_v(fFieldComponents.z()));
  }

  /** @brief Field value setter */
  void SetFieldValue(const Vector3D<float>& fieldValue) { fFieldComponents = fieldValue; }

  /** @brief Field value getter */
  vecgeom::Vector3D<float> GetConstantFieldValue() const { return fFieldComponents; }

private:
  vecgeom::Vector3D<float> fFieldComponents;
};

UniformMagField::UniformMagField(double vField,
                                 double vTheta,
                                 double vPhi)
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
#endif
