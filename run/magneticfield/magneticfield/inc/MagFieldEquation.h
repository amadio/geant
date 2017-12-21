// Approach is derived from the Geant4 class G4MagFieldEquation
//

#include <cmath>
#include <iostream>

#ifndef MagFieldEquation_H
#define MagFieldEquation_H  1

#include "base/Vector3D.h"
#include <Geant/Config.h>
#include <Geant/VectorTypes.h>

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

// #define OUTSIDE_MagFieldEquation 1

template <class Field>
class MagFieldEquation
{
public:
//  static const unsigned int  N   = Size;
  using Double_v = Geant::Double_v;
  using Float_v = Geant::Float_v;
  
  static constexpr double gCof   = geant::kCLight;  //   / fieldUnits::meter ;

  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  // Expected constant value:
  // static constexpr double gCof    = Constants::c_light * fieldUnits::second /
  //     ( 1.0e9 * fieldUnits::meter * fieldUnits::meter );

  MagFieldEquation(Field* pF) : fPtrField(pF) {}

  MagFieldEquation(const MagFieldEquation& right) : fPtrField(right.fPtrField) {}

  ~MagFieldEquation() {}

  GEANT_FORCE_INLINE
  Field* GetField() const { return fPtrField; } 
  
  template <typename Real_v>
  GEANT_FORCE_INLINE
  void RightHandSide( const Real_v y[],
                            Real_v charge,
                            Real_v dydx[] ) const
#ifdef OUTSIDE_MagFieldEquation
     ;
#else  
  {
    Vector3D<Real_v> Bfield;
    FieldFromY( y, Bfield );
    EvaluateRhsGivenB( y, Bfield, charge, dydx );
  }
#endif  
  template <typename Real_v>
  GEANT_FORCE_INLINE
  void EvaluateRhsGivenB(const Real_v            y[],
                         const Vector3D<Real_v> &B,
                         const Real_v           &charge,
                               Real_v            dydx[] ) const
  {
    // ThreeVectorD momentum( y[3], y[4], y[5]);
    Real_v momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    Real_v inv_momentum_magnitude = Real_v(1.) / vecCore::math::Sqrt( momentum_mag_square );
    // Real_v inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);
  
    Real_v cof = charge * Real_v(gCof) * inv_momentum_magnitude;
    
    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)
  }

  template <typename Real_v>
  GEANT_FORCE_INLINE
  void FieldFromY(const Real_v y[], Vector3D<Real_v> &Bfield ) const
  {
    fPtrField->GetFieldValue( Vector3D<Real_v>(y[0], y[1], y[2]), Bfield );  
  }

  template <typename Real_v>
  void PrintInputFieldAndDyDx(const Real_v y[], Real_v charge, Real_v dydx[] ) const
  {
    RightHandSide(y, dydx);

    // Obtain the field value
    Vector3D<Real_v>  Bfield;
    FieldFromY( y, Bfield );
    EvaluateRhsGivenB(y, charge, Bfield, dydx);

    std::cout.precision(8);

    std::cout << "\n# Input & B field \n";
    std::cout.setf (std::ios_base::scientific);
    std::cout << " Position = " << y[0] << " " << y[1] << " " << y[2] << std::endl;
    std::cout << " Momentum = " << y[3] << " " << y[4] << " " << y[5] << std::endl;
    std::cout << " B-field  = " << Bfield[0] << " " << Bfield[1] << " " << Bfield[2] << std::endl;
    std::cout.unsetf(std::ios_base::scientific);

    std::cout << "\n# 'Force' from B field \n";
    std::cout.setf (std::ios_base::fixed);
    std::cout << " dy/dx [0-2] (=dX/ds) = " << dydx[0]   << " " << dydx[1]   << " " << dydx[2] << std::endl;
    std::cout << " dy/dx [3-5] (=dP/ds) = " << dydx[3]   << " " << dydx[4]   << " " << dydx[5] << std::endl;
    std::cout.unsetf(std::ios_base::fixed);
  }

  template <typename Real_v>
  void PrintAll( Real_v const  y[],
                 const Vector3D<Real_v> &B,
                 Real_v        charge,
                 Real_v        cof,                     
                 Real_v const  dydx[]  ) const
  {
    using ThreeVector = Vector3D<Real_v>;
    using geant::kilogauss;
    
    std::cout.precision(8);
    std::cout << "Equation:  gCof= " << gCof << " charge= " << charge << " cof= " << cof
              << " Bfield= " << B << std::endl;
    std::cout << "            dx/ds  = " << dydx[0] << " " << dydx[1] << " " << dydx[2]
              << " - mag= " << std::sqrt( dydx[0] * dydx[0] + dydx[1] * dydx[1] + dydx[2] * dydx[2] ) << std::endl;
    std::cout << "            dp/ds  = " << dydx[3] << " " << dydx[4] << " " << dydx[5]
              << " - mag= " << std::sqrt( dydx[3] * dydx[3] + dydx[4] * dydx[4] + dydx[5] * dydx[5] ) << std::endl;

    Real_v Bmag = ThreeVector(B[0], B[1], B[2]).Mag();
    std::cout << "            B-field= " << B[0] / kilogauss << " " << B[1] / kilogauss << " "
              << B[2] / kilogauss << "  ( KGaus ) mag= " << Bmag << std::endl;
    std::cout << "               P  = " << y[3] << " " << y[4] << " " << y[5]
              << " = mag= " << ThreeVectorD(y[3],y[4],y[5]).Mag() << std::endl;
  }
  

private:
  enum { G4maximum_number_of_field_components = 24 };
  Field    *fPtrField = nullptr; // The field object
};

#ifdef OUTSIDE_MagFieldEquation
template <typename Real_v>
GEANT_FORCE_INLINE
void
template <class Field>
  MagFieldEquation<Field>::
       RightHandSide( const Real_v y[],
                            Real_v charge,
                            Real_v dydx[] ) const
  {
    Vector3D<Real_v> Bfield;
    FieldFromY( y, Bfield );
    EvaluateRhsGivenB( y, Bfield, charge, dydx );
  }
#endif

#undef OUTSIDE_MagFieldEquation

#endif  // MagFieldEquation_H
