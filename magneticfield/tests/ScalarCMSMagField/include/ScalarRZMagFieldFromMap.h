#ifndef  SCALARRZMAGFIELDFROMMAP_H_
#define  SCALARRZMAGFIELDFROMMAP_H_

#include <vector>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"

#include <Geant/SystemOfUnits.h>

#include "Geant/VVectorField.h"

// typedef vecgeom::Vector3D<double> Vector3D;
// typedef vecgeom::SOA3D<double> SOA3D;

class ScalarRZMagFieldFromMap : public VVectorField {
public:
  ScalarRZMagFieldFromMap();
  ScalarRZMagFieldFromMap(std::string inputMap);
  ScalarRZMagFieldFromMap(const ScalarRZMagFieldFromMap &right);

  // New stuff
  // Takes as input x,y,z; Gives output Bx,By,Bz
  void GetFieldValueXYZ(const vecgeom::Vector3D<double> &position, vecgeom::Vector3D<double> &xyzField);

  // Stores rz field as well for cross checking purpose
  void GetFieldValueTest(const vecgeom::Vector3D<double> &position, vecgeom::Vector3D<double> &rzField);

  // @brief   evaluate field value at input position (SOA3D)
  void ObtainFieldValue(const vecgeom::Vector3D<double> &position,
                        vecgeom::Vector3D<double> &field) override final
  { GetFieldValueXYZ(position, field); }

  virtual void ObtainFieldValueSIMD(const Vector3D<Double_v> &position, Vector3D<Double_v> &fieldValue) override final;

  // @brief   extra method for an arbitrary length SOA3D of positions
  void GetFieldValues(const vecgeom::SOA3D<double> &posVec, vecgeom::SOA3D<double> &fieldVec);
  
  // Reads data from given 2D magnetic field map. Can be easily modified to read a given 2D map, in case the file
  // changes
  void ReadVectorData(std::string inputMap);

  ~ScalarRZMagFieldFromMap();

  /** @brief For old interface - when cloning was needed for each thread */
  ScalarRZMagFieldFromMap* Copy() const { return new ScalarRZMagFieldFromMap(*this); }  
   
  virtual VVectorField* Clone() const override final { return Copy(); }

  ScalarRZMagFieldFromMap* CloneOrSafeSelf(bool *isSafe) {
     if( isSafe ) *isSafe = false;
     return Copy();
  }

  
private:
  //  Invariants -- parameters of the field
  const double millimeter = 1.0; // Currently -- to be native GeantV

  const double kRMax     = 9000. * millimeter;  //  Maximum value of R =  9.00 meters
  const double kZMax     = 16000. * millimeter; //  Max value of Z = 16.00 meters
  const int kNoZValues   = 161;
  const int kNoRValues   = 181;
  const int kHalfZValues = 80;

  // Derived values
  // kRDiff and kZDiff take care of mm because they come from kRMax and kZMax which have mm in them
  const double kRDiff = kRMax / (kNoRValues - 1);     //  Radius increment between lattice points
  const double kZDiff = 2 * kZMax / (kNoZValues - 1); //  Z increment

  const double kZ0       = -kZMax;
  const double kRDiffInv = 1.0 / kRDiff;
  const double kZDiffInv = 1.0 / kZDiff;
  const double kAInverse = geant::units::tesla / (kRDiff * kZDiff);  // Values in file assumed in tesla

  // For (R,Z) pairs : gives field in cylindrical coordinates in rzfield
  void GetFieldValueRZ(const double radius, const double z, vecgeom::Vector3D<double> &rzField);
  void GetFieldValueRZ(std::vector<double> radius, std::vector<double> z);

  // Used to convert cartesian coordinates to cylindrical coordinates R-Z-phi
  //  Does not calculate phi
  inline void CartesianToCylindrical(const vecgeom::Vector3D<double> &cart, double cyl[2]);

  // Converts cylindrical magnetic field to field in cartesian coordinates
  inline void CylindricalToCartesian(const vecgeom::Vector3D<double> &B_cyl, const double sinTheta, const double cosTheta,
                                     vecgeom::Vector3D<double> &B_cart);

private:
  std::vector<double> fRadius, fPhi, fZ, fBr, fBz, fBphi;
};

inline void ScalarRZMagFieldFromMap::CartesianToCylindrical(const vecgeom::Vector3D<double> &cart, double cyl[2])
{
  // cyl[3] =[r,z,phi]
  cyl[0] = sqrt(cart[0] * cart[0] + cart[1] * cart[1]); // r = sqrt(x^2 + y^2)
  cyl[1] = cart[2];                                     // z = z
}

inline void ScalarRZMagFieldFromMap::CylindricalToCartesian(const vecgeom::Vector3D<double> &rzField, const double sinTheta, const double cosTheta,
                                             vecgeom::Vector3D<double> &xyzField)
{
  // B_cyl[] has r, phi and z
  xyzField[0] = rzField[0] * cosTheta - rzField[1] * sinTheta; // Bx= Br cos(theta) - Bphi sin(theta)
  xyzField[1] = rzField[0] * sinTheta + rzField[1] * cosTheta; // By = Br sin(theta) + Bphi cos(theta)
  xyzField[2] = rzField[2];                                    // Bz = Bz
}

#endif
