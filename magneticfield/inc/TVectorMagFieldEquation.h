// Approach is derived from the Geant4 class G4MagFieldEquation
//

#include <cmath>

#include "GUVEquationOfMotion.h"
#include "GUVVectorEquationOfMotion.h"
#include "TMagFieldEquation.h"

//  Ensure that equation Right Hand Side is inlined - may be compiler dependend
#define INLINERHS 1

#ifdef  INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline)) 
#else
#define REALLY_INLINE   inline
#endif

#ifndef TVECTORMAGFIELDEQUATION_H
#define TVECTORMAGFIELDEQUATION_H  1

// #include <vector>
#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "Units.h"
#include "Constants.h"
//  Update to GeantV units ASAP

template <class Field, unsigned int Size>
class TVectorMagFieldEquation :  public GUVVectorEquationOfMotion
{
   using Double_v = Geant::Double_v;
   //using Float_v  = Geant::Float_v;

   template <typename T>
   using Vector3D = vecgeom::Vector3D<T>;
   
   public:
//     static const unsigned int  N   = Size;
     // static const double fCof   = Constants::c_light;   // Was constexpr

     TVectorMagFieldEquation(Field* pF) : GUVVectorEquationOfMotion(pF) { fPtrField = pF; }

     TVectorMagFieldEquation(const TVectorMagFieldEquation& );
     ~TVectorMagFieldEquation()  {}  // Was virtual - but now no inheritance

     TVectorMagFieldEquation<Field,Size>* Clone() const;
     TVectorMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool& safe);
     TVectorMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool* safe=0);
     
     REALLY_INLINE  
     void GetFieldValue(const Double_v Point[4],
                              Double_v Value[]) const
     {
        fPtrField->GetFieldValue(Point, Value);
     }

     inline 
     void RightHandSide(const Double_v y[], 
                        const Double_v charge, 
                              Double_v dydx[] ) const;

     REALLY_INLINE
     void TEvaluateRhsGivenB( const Double_v y[],
                              const Vector3D<Double_v> B,  // Was double B[3],
                              const Double_v charge,
                                    Double_v dydx[] ) const;

     // virtual
     void EvaluateRhsGivenB( const Double_v y[],
                             const Vector3D<Double_v> B,  // Was const double B[3],
                             const Double_v charge,
                                   Double_v dydx[] ) const override final
     { TEvaluateRhsGivenB( y, B, charge, dydx); }

     REALLY_INLINE
     void FieldFromY(const Double_v y[], 
                           Double_v Bfield[] ) const;

     REALLY_INLINE
     void FieldFromY(const Double_v y[], 
                           Vector3D<Double_v> &Bfield ) const;

     REALLY_INLINE
     void PrintInputFieldAndDyDx(const Double_v y[],  
                                 const Double_v charge,  
                                       Double_v dydx[] ) const;

   private:
     enum { G4maximum_number_of_field_components = 24 };
     Field *fPtrField;
     double   fParticleCharge;
};

template <class Field, unsigned int Size>
   TVectorMagFieldEquation<Field,Size>
   ::TVectorMagFieldEquation(const TVectorMagFieldEquation& right)
   :  GUVVectorEquationOfMotion( (GUVField*) 0 ),
      fPtrField( right.fPtrField->CloneOrSafeSelf( (bool *)0 ) )
      // fPtrField( new Field(right.fPtrField) )
{
   // G4bool threadSafe;
   // fPtrField = right.fPtrField->CloneOrSafeSelf( &threadSafe );

   // std::cout <<  "TVectorMagFieldEquation - copy constructor called." << std::endl;
   GUVVectorEquationOfMotion::SetFieldObj( fPtrField ); //  Also stored in base class ... for now
}

template <class Field, unsigned int Size>
   TVectorMagFieldEquation<Field,Size>*
   TVectorMagFieldEquation<Field,Size>
   ::CloneOrSafeSelf(bool& safe)
{
   TVectorMagFieldEquation<Field,Size>* equation;
   Field* pField=
      fPtrField->CloneOrSafeSelf(safe);

   std::cerr << " #TVectorMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool& safe) called# " << std::endl;

      equation = new TVectorMagFieldEquation( pField );
      safe= false;

}



template <class Field, unsigned int Size>
REALLY_INLINE
   void  TVectorMagFieldEquation<Field, Size>
   ::TEvaluateRhsGivenB( const Double_v             y[],
                         const Vector3D<Double_v>   B,
                         const Double_v             charge,
                         Double_v                   dydx[]  ) const
{
  
    Double_v momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    Double_v inv_momentum_magnitude = 1. / vecCore::math::Sqrt( momentum_mag_square );

    Double_v cof = charge * ((Double_v) Constants::c_light)   // Was fCof
                   * inv_momentum_magnitude;

    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

}

template <class Field, unsigned int Size>
REALLY_INLINE
void
TVectorMagFieldEquation<Field,Size>
   ::FieldFromY(const Double_v y[], 
                      Double_v Bfield[3] ) const
{
    // double  Bfield[3];  //G4maximum_number_of_field_components];
    Double_v PositionAndTime[4];
    PositionAndTime[0] = y[0];
    PositionAndTime[1] = y[1];
    PositionAndTime[2] = y[2];
    PositionAndTime[3] = 0;        // Time

    GetFieldValue(PositionAndTime, Bfield) ;
}

template <class Field, unsigned int Size>
REALLY_INLINE
void
TVectorMagFieldEquation<Field,Size>
   ::FieldFromY(const Double_v             y[],  
                      Vector3D<Double_v>   &Bfield ) const
{
    Vector3D<Double_v> Position( y[0], y[1], y[2] );

    fPtrField->Field::GetFieldValue( Position, Bfield );
}


template <class Field, unsigned int Size>
REALLY_INLINE
void
TVectorMagFieldEquation<Field,Size>
   ::RightHandSide(const Double_v y[], 
                   const Double_v charge, 
                         Double_v dydx[] ) const
{
    Vector3D<Double_v> BfieldVec;

    FieldFromY( y, BfieldVec );
    TEvaluateRhsGivenB( y, BfieldVec, charge, dydx );
}

#include <iostream>   // For debuging only
using std::cout;
using std::endl;

template <class Field, unsigned int Size>
REALLY_INLINE
void
TVectorMagFieldEquation<Field,Size>
   ::PrintInputFieldAndDyDx(const Double_v y[], 
                            const Double_v charge, 
                                  Double_v dydx[] ) const
{

    RightHandSide(y, dydx);

    // Obtain the field value
    Double_v  Bfield[3];  //G4maximum_number_of_field_components];
    FieldFromY( y, charge, Bfield );
    TEvaluateRhsGivenB(y, Bfield, charge, dydx);

    cout.precision(8);

    // cout.setf (std::ios_base::fixed);
    // cout << " Position = " << y[0] << " " << y[1] << " " << y[3] << endl;
    // cout.unsetf(std::ios_base::fixed);
    cout << "\n# Input & B field \n";
    cout.setf (std::ios_base::scientific);
    cout << " Position = " << y[0] << " " << y[1] << " " << y[2] << endl;
    cout << " Momentum = " << y[3] << " " << y[4] << " " << y[5] << endl;
    cout << " B-field  = " << Bfield[0] << " " << Bfield[1] << " " << Bfield[2] << endl;
    cout.unsetf(std::ios_base::scientific);

    cout << "\n# 'Force' from B field \n";
    cout.setf (std::ios_base::fixed);
    cout << " dy/dx [0-2] (=dX/ds) = " << dydx[0]   << " " << dydx[1]   << " " << dydx[2] << endl;
    cout << " dy/dx [3-5] (=dP/ds) = " << dydx[3]   << " " << dydx[4]   << " " << dydx[5] << endl;
    cout.unsetf(std::ios_base::fixed);
}
#endif  // TMAGFIELDEQUATION_H
