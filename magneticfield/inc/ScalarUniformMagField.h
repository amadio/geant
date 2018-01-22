//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef ScalarUniformMagField_H
#define ScalarUniformMagField_H

#include "VVectorField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "SystemOfUnits.h"
// #include "PhysicalConstants.h"

class ScalarUniformMagField : public VVectorField
{
public:

  static constexpr int   gNumFieldComponents= 3;
  static constexpr bool  gFieldChangesEnergy= false;
  
  /** @brief Constructor providing the constant field value (cartesian) */
  ScalarUniformMagField( const vecgeom::Vector3D<float>& fieldVector )
    : VVectorField(gNumFieldComponents, gFieldChangesEnergy),
      fFieldComponents(fieldVector) {}

  /** @brief Constructor providing the constant field value (spherical) */
  ScalarUniformMagField(double vField,
                   double vTheta,
                   double vPhi  );

  /** @brief Destructor */
  ~ScalarUniformMagField() {}

  /** @brief Copy constructor */
  ScalarUniformMagField(const ScalarUniformMagField &p)
    : VVectorField(gNumFieldComponents, gFieldChangesEnergy),
      fFieldComponents(p.fFieldComponents) {}

  /** Assignment operator */
  ScalarUniformMagField& operator = (const ScalarUniformMagField &p);

  /** @brief Scalar interface for field retrieval */
  virtual void  ObtainFieldValue( const Vector3D<double> &position,
                                     Vector3D<double>  &fieldValue )
  {
    GetFieldValue<double>(position, fieldValue);
  }

  /** @brief Vector interface for field retrieval */
  virtual void  ObtainFieldValueSIMD( const Vector3D<Double_v> &position, 
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

  /** IS THIS NEEDED ? */
  ScalarUniformMagField* Clone() const {  return new ScalarUniformMagField( *this ); }

  ScalarUniformMagField* CloneOrSafeSelf( bool* pSafe )
  {
     if( pSafe ) *pSafe= true;
     return this; // ->CloneOrSafeSelf(*pSafe);
  }
  //  Class is thread-safe, can use 'self' instead of clone

private:
  vecgeom::Vector3D<float> fFieldComponents;
};

ScalarUniformMagField& 
ScalarUniformMagField::operator = (const ScalarUniformMagField &p)
   // Copy constructor and assignment operator.
{
   if (&p == this) return *this;
   // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
   fFieldComponents = p.fFieldComponents;
   return *this;
}

ScalarUniformMagField::ScalarUniformMagField(double vField,
                                   double vTheta,
                                   double vPhi     )
  : VVectorField(gNumFieldComponents, gFieldChangesEnergy)
{
   if ( (vField<0) || (vTheta<0) || (vTheta>geant::kPi) || (vPhi<0) || (vPhi>geant::kTwoPi) )
   {
      // Exception("ScalarUniformMagField::ScalarUniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in ScalarUniformMagField::ScalarUniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << geant::kPi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << geant::kTwoPi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << geant::kTwoPi << std::endl;
   }
   fFieldComponents.Set( vField*std::sin(vTheta)*std::cos(vPhi),
                         vField*std::sin(vTheta)*std::sin(vPhi),
                         vField*std::cos(vTheta)                );
}
#endif
