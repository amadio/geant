// Approach is derived from the Geant4 class G4MagFieldEquation
//

#include <cmath>

#include "VScalarEquationOfMotion.h"

#ifndef SCALAR_MAG_FIELD_EQUATION_H
#define SCALAR_MAG_FIELD_EQUATION_H  1

// #include <vector>
#include "base/Vector3D.h"

#include "Units.h"
#include "Constants.h"
//  Update to GeantV units ASAP

template <class Field, unsigned int Size>
class ScalarMagFieldEquation : public VScalarEquationOfMotion
{
public:
//  static const unsigned int  N   = Size;
  static constexpr double gCof   = Constants::c_light;  //   / fieldUnits::meter ;

  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  // Expected constant value:
  // static constexpr double gCof    = Constants::c_light * fieldUnits::second /
  //     ( 1.0e9 * fieldUnits::meter * fieldUnits::meter );

  ScalarMagFieldEquation(Field* pF) : VScalarEquationOfMotion(pF) { fPtrField = pF; }

  ScalarMagFieldEquation(const ScalarMagFieldEquation& );
  ~ScalarMagFieldEquation()  {}  // Was virtual - but now no inheritance

  ScalarMagFieldEquation<Field,Size>* Clone() const;
  ScalarMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool& safe);
  ScalarMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool* safe=0);
    // If class is thread safe, return self ptr.  Else return clone

  Field* GetField() { return fPtrField; } 

  inline // GEANT_FORCE_INLINE
  void RightHandSide(const          double   y[],
                     const Vector3D<double> &position,
                                    double   charge,
                                    double   dydx[],
                           Vector3D<double> &BfieldVec) const;

  inline // GEANT_FORCE_INLINE
  void RightHandSide(const          double   y[],
                                    double   charge,
                                    double   dydx[],
                           Vector3D<double> &BfieldVec) const;

  inline // GEANT_FORCE_INLINE
  void RightHandSide(const          double   y[],
                                    double   charge,
                                    double   dydx[]) const;

  virtual
  void EvaluateRhsGivenB( const          double   y[],
                          const Vector3D<double> &B,
                                         double   charge,
                                         double   dydx[] ) const;
  
  template <typename Real_v>
  void EvaluateRhsGivenB(const Real_v            y[],
                         const Vector3D<Real_v> &B,
                         const Real_v           &charge,
                               Real_v            dydx[] ) const
  {
    // ThreeVectorD momentum( y[3], y[4], y[5]);
    double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    double inv_momentum_magnitude = 1. / std::sqrt( momentum_mag_square );
    // double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);
  
    double cof = charge*gCof*inv_momentum_magnitude;
    
    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)
  }

  GEANT_FORCE_INLINE
  void TEvaluateRhsReturnB( const double    y[],
                            double          dydx[],
                            double          charge,
                            vecgeom::Vector3D<double> &field  ) const;
  
  GEANT_FORCE_INLINE
  void FieldFromY(const double y[], vecgeom::Vector3D<double> &Bfield ) const;

  GEANT_FORCE_INLINE
  void PrintInputFieldAndDyDx(const double y[], double charge, double dydx[] ) const;

  GEANT_FORCE_INLINE
  void PrintAll(  double const  y[],
                 const vecgeom::Vector3D<double> &B,
                 double        charge,
                 double        cof,                     
                 double const  dydx[]  ) const;

private:
  enum { G4maximum_number_of_field_components = 24 };
  Field    *fPtrField;
};

template <class Field, unsigned int Size>
   ScalarMagFieldEquation<Field,Size>::ScalarMagFieldEquation(const ScalarMagFieldEquation& right)
     : VScalarEquationOfMotion( (VVectorField*) 0 ),
       fPtrField( right.fPtrField->CloneOrSafeSelf( (bool *)0 ) )
      // fPtrField( new Field(right.fPtrField) )
{
   // G4bool threadSafe;
   // fPtrField = right.fPtrField->CloneOrSafeSelf( &threadSafe );

   // std::cout <<  "ScalarMagFieldEquation - copy constructor called." << std::endl;
   VScalarEquationOfMotion::SetFieldObj( fPtrField ); //  Also stored in base class ... for now
}

template <class Field, unsigned int Size>
   ScalarMagFieldEquation<Field,Size> *ScalarMagFieldEquation<Field,Size>::Clone() const
{
   // bool safe= false;  // Field* pField= fPtrField->CloneOrSafeSelf(safe);
   Field* cloneField= fPtrField->Clone();   
   std::cerr << " #ScalarMagFieldEquation<Field,Size>::Clone() called# " << std::endl;
   return new ScalarMagFieldEquation( cloneField );
}

template <class Field, unsigned int Size>
ScalarMagFieldEquation<Field,Size> *ScalarMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool& safe)
{
   ScalarMagFieldEquation<Field,Size>* equation;
   Field* pField=
      fPtrField->CloneOrSafeSelf(safe);
   // If Field does not have such a method:
   //  = new Field( fPtrField ); // Need copy constructor.
   //  safe= false;

   std::cerr << " #ScalarMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool& safe) called# " << std::endl;

   // safe = safe && fClassSafe;
   // if( safe )  {  equation = this; }
   //    Can be introduced when Equation is thread safe -- no state
   //     --> For now the particle Charge is preventing this 23.11.2015
   // else {
      equation = new ScalarMagFieldEquation( pField );
      safe= false;
   // }

   return equation;   
}

template <class Field, unsigned int Size>
ScalarMagFieldEquation<Field,Size> *ScalarMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool* pSafe)
{
   bool safeLocal;
   std::cerr << " #ScalarMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool* safe) called#" << std::endl;
   if( !pSafe ) pSafe= &safeLocal;
   auto equation= CloneOrSafeSelf( pSafe );
   return equation;   
}


template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void  ScalarMagFieldEquation<Field, Size>::EvaluateRhsGivenB(
                          const          double   y[],
                          const Vector3D<double> &B,
                                         double   charge,
                                         double   dydx[]  ) const
{
    // ThreeVectorD momentum( y[3], y[4], y[5]);
    double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    double inv_momentum_magnitude = 1. / std::sqrt( momentum_mag_square );
    // double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);
  
    double cof = charge*gCof*inv_momentum_magnitude;
    
    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

}

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void  ScalarMagFieldEquation<Field, Size>::PrintAll(
                double const  y[],
                const vecgeom::Vector3D<double> &B,
                double        charge,
                double        cof,
                double const  dydx[]  ) const
{
    using ThreeVectorD = vecgeom::Vector3D<double>;
   
    printf("Equation:  gCof = %8.4g  charge= %f cof= %10.5g   B-field= %f %f %f \n",
           gCof, charge, cof, B[0], B[1], B[2] );
    // printf("               X  = %12.6g %12.6g %12.6g - mag %12.6g\n",  y[0], y[1], y[2] );

    printf("            dx/ds  = %12.6g %12.6g %12.6g - mag %12.6g\n",
           dydx[0], dydx[1], dydx[2],  std::sqrt( dydx[0] * dydx[0] + dydx[1] * dydx[1] + dydx[2] * dydx[2] ) );
    printf("            dp/ds  = %12.6g %12.6g %12.6g - mag %12.6g\n",
           dydx[3], dydx[4], dydx[5],  std::sqrt( dydx[3] * dydx[3] + dydx[4] * dydx[4] + dydx[5] * dydx[5] ) );

    double Bmag2chk = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
    printf("            B-field= %10.3f %10.3f %10.3f  ( KGaus ) mag= %10.4f\n",
           B[0] / fieldUnits::kilogauss ,
           B[1] / fieldUnits::kilogauss ,
           B[2] / fieldUnits::kilogauss ,
           std::sqrt(Bmag2chk) );

    printf("               P  = %12.6g %12.6g %12.6g - mag %12.6g\n",  y[3], y[4], y[5],
           ThreeVectorD(y[3],y[4],y[5]).Mag() );

    return ;
}

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>
  ::FieldFromY(const double y[],  vecgeom::Vector3D<double> &Bfield ) const
{
   fPtrField->GetFieldValue( Vector3D<double>(y[0], y[1], y[2]), Bfield );
}


template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>
  ::RightHandSide(const double y[], double charge, double dydx[] ) const
{

    vecgeom::Vector3D<double> BfieldVec;

    FieldFromY( y, BfieldVec );
    EvaluateRhsGivenB( y, BfieldVec, charge, dydx );
}

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>::RightHandSide(
                    const double  y[],
                    double  charge,
                    double  dydx[],
                    Vector3D<double> &BfieldVec ) const
{
    FieldFromY( y, BfieldVec );
    EvaluateRhsGivenB( y, BfieldVec, charge, dydx );
}

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>
  ::RightHandSide(const                   double   y[],
                  const          Vector3D<double> &position,
                                          double   charge,
                                          double   dydx[],
                                 Vector3D<double> &BfieldVec) const
{
   fPtrField->GetFieldValue( position, BfieldVec );
   EvaluateRhsGivenB( y, BfieldVec, charge, dydx );
}

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>
  ::TEvaluateRhsReturnB( const double           y[],
                          double          dydx[],
                          double          charge,
                    vecgeom::Vector3D<double> &field
                  ) const
{
   using ThreeVector = vecgeom::Vector3D<double>;
   
   GetFieldValue( ThreeVector(y[0], y[1], y[2]), field) ;
   EvaluateRhsGivenB( y, field, charge, dydx );
}

#include <iostream>   // For debuging only

template <class Field, unsigned int Size>
GEANT_FORCE_INLINE
void ScalarMagFieldEquation<Field,Size>
  ::PrintInputFieldAndDyDx(const double y[], double charge, double dydx[] ) const
{
  RightHandSide(y, dydx);

  // Obtain the field value
  Vector3D<double>  Bfield;
  FieldFromY( y, Bfield );
  EvaluateRhsGivenB(y, charge, Bfield, dydx);

  std::cout.precision(8);

  // cout.setf (std::ios_base::fixed);
  // cout << " Position = " << y[0] << " " << y[1] << " " << y[3] << std::endl;
  // cout.unsetf(std::ios_base::fixed);
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
#endif  // SCALAR_MAG_FIELD_EQUATION_H
