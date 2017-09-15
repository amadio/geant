//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TUniformMagField_H
#define TUniformMagField_H

#include "GUVField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

class TUniformMagField : public GUVField
{
public:

  static constexpr int   gNumFieldComponents= 3;
  static constexpr bool  gFieldChangesEnergy= false;
  
  /** @brief Constructor providing the constant field value (cartesian) */
  TUniformMagField( const vecgeom::Vector3D<float>& fieldVector )
    : GUVField(gNumFieldComponents, gFieldChangesEnergy),
      fFieldComponents(fieldVector) {}

  /** @brief Constructor providing the constant field value (spherical) */
  TUniformMagField(double vField,
                   double vTheta,
                   double vPhi  );

  /** @brief Destructor */
  ~TUniformMagField() {}

  /** @brief Copy constructor */
  TUniformMagField(const TUniformMagField &p)
    : GUVField(gNumFieldComponents, gFieldChangesEnergy),
      fFieldComponents(p.fFieldComponents) {}

  /** Assignment operator */
  TUniformMagField& operator = (const TUniformMagField &p);

  /** @brief Scalar interface for field retrieval */
  virtual void  GetFieldValue( const Vector3D<double> &position, 
                                     Vector3D<double>  &fieldValue )
  {
    GetFieldValue<double>(position, fieldValue);
  }

  /** @brief Vector interface for field retrieval */
  virtual void GetFieldValueSIMD( const Vector3D<Double_v> &position, 
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
  TUniformMagField* Clone() const {  return new TUniformMagField( *this ); }

  TUniformMagField* CloneOrSafeSelf( bool* pSafe )
  {
     if( pSafe ) *pSafe= true;
     return this; // ->CloneOrSafeSelf(*pSafe);
  }
  //  Class is thread-safe, can use 'self' instead of clone

private:
  vecgeom::Vector3D<float> fFieldComponents;
};

TUniformMagField& 
TUniformMagField::operator = (const TUniformMagField &p)
   // Copy constructor and assignment operator.
{
   if (&p == this) return *this;
   // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
   fFieldComponents = p.fFieldComponents;
   return *this;
}

TUniformMagField::TUniformMagField(double vField,
                                   double vTheta,
                                   double vPhi     )
  : GUVField(gNumFieldComponents, gFieldChangesEnergy)
{
   if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
   {
      // Exception("TUniformMagField::TUniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in TUniformMagField::TUniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << Constants::twopi << std::endl;
   }
   fFieldComponents.Set( vField*std::sin(vTheta)*std::cos(vPhi),
                         vField*std::sin(vTheta)*std::sin(vPhi),
                         vField*std::cos(vTheta)                );
}
#endif
