//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef UniformMagField_H
#define UniformMagField_H

#include <iostream>

#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "VVectorField.h"

class UniformMagField : public VVectorField
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
     :  VVectorField( 3, true ),  // Field does not change energy
        fFieldComponents(fieldVector) {}

  /** @brief Constructor providing the constant field value (spherical) */
  UniformMagField( double vField,
                   double vTheta,
                   double vPhi  );

  /** @brief Destructor */
  ~UniformMagField() {}

  /** @brief Templated field interface */
  template <typename Real_v>
  void GetFieldValue( const Vector3D<Real_v> & /*position*/,
                            Vector3D<Real_v> &fieldValue )
  {
    fieldValue.Set(Real_v(fFieldComponents.x()),
                   Real_v(fFieldComponents.y()),
                   Real_v(fFieldComponents.z()));
  }
  
  /** @brief Fast Scalar interface for field retrieval */
  void  GetFieldValue( const Vector3D<double> &position, 
                             Vector3D<double> &fieldValue )
  {
    GetFieldValue<double>(position, fieldValue);
  }

  /** @brief Fast Vector interface for field retrieval */
  void GetFieldValueSIMD( const Vector3D<Double_v> &position, 
                                Vector3D<Double_v> &fieldValue )
  {
    GetFieldValue<Double_v>(position, fieldValue);
  }

  /** @brief Scalar interface for field retrieval  */
  virtual void ObtainFieldValue( const Vector3D<double> &position,
                                       Vector3D<double>  &fieldValue );

  /** @brief Vector interface for field retrieval */
  virtual void ObtainFieldValueSIMD( const Vector3D<Double_v> &position, 
                                           Vector3D<Double_v> &fieldValue );
  
  /** @brief Field value setter */
  void SetFieldValue(const Vector3D<float>& fieldValue) { fFieldComponents = fieldValue; }

  /** @brief Field value getter */
  vecgeom::Vector3D<float> GetConstantFieldValue() const { return fFieldComponents; }

  // STATE
private:
  vecgeom::Vector3D<float> fFieldComponents;
};


#endif
