// Approach is derived from the Geant4 class G4MagFieldEquation
// 

#include <cmath>

#include "GUVEquationOfMotion.h"

//  Ensure that equation Right Hand Side is inlined - may be compiler dependend
#define INLINERHS 1

#ifdef  INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline)) 
#else
#define REALLY_INLINE   inline
#endif

#ifndef TMAGFIELDEQUATION_H
#define TMAGFIELDEQUATION_H  1

template 
<class Field, unsigned int Size>
class TMagFieldEquation : public GUVEquationOfMotion
{
   public:
     typedef Field T_Field;
     static const unsigned int  N   = Size;
     static constexpr double meter   = 100.0; // ASSUME centimeter native unit
     static constexpr double second  = 1.0;   // ASSUME second  is native unit
     static constexpr double c_light = 2.99792458e+8 * meter/second; // Units TBC
     static constexpr double eplus   = 1.0 ;   // Units TBC
     static constexpr double fCof    = eplus * c_light ;
   
     TMagFieldEquation(T_Field* pF) : GUVEquationOfMotion(pF) { fPtrField = pF; }
     ~TMagFieldEquation()  {}  // Was virtual - but now no inheritance

     REALLY_INLINE  // inline __attribute__((always_inline))     
     void GetFieldValue(const double Point[4],
                              double Value[]) const
     {
        fPtrField->T_Field::GetFieldValue(Point, Value);
     }

     REALLY_INLINE
        void RightHandSide(const double y[], /*double charge,*/ double dydx[] ) const;

     REALLY_INLINE
     void TEvaluateRhsGivenB( const double y[],
                              const double B[3],
                                //  double charge, 
                             double dydx[] ) const;

     // virtual
     void EvaluateRhsGivenB( const double y[],
                             const double B[3],
                              //   double charge, 
                             double dydx[] ) const
     { TEvaluateRhsGivenB( y, B, /*charge,*/ dydx); }

     void InitializeCharge(double particleCharge) final
      { fParticleCharge= particleCharge;  GUVEquationOfMotion::InformReady();  }
      
      void InvalidateParameters() final { GUVEquationOfMotion::InformDone();}

   private:
     enum { G4maximum_number_of_field_components = 24 };
     T_Field *fPtrField;
     double    fParticleCharge;  // Can be moved to concrete Equation (as a type!)
                                 // - to generalised for other forces
};

template 
<class Field, unsigned int Size>
REALLY_INLINE
   void  TMagFieldEquation<Field, Size>
   ::TEvaluateRhsGivenB( const double y[],
                         const double B[3],
                           //  double charge,
                               double dydx[]  ) const
{
    double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    double inv_momentum_magnitude = 1. / std::sqrt( momentum_mag_square );
    // double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);

    double cof = fParticleCharge*fCof*inv_momentum_magnitude;

    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

    return ;
}

template 
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size> 
   ::RightHandSide(const double y[], /* double charge, */ double dydx[] ) const
{
    double Point[4];  //G4maximum_number_of_field_components]; 
    double  PositionAndTime[4];
    PositionAndTime[0] = y[0];
    PositionAndTime[1] = y[1];
    PositionAndTime[2] = y[2];
    PositionAndTime[3] = 0;        // Time
    // PositionAndTime[3] = y[7];  //  --> extersion?
    GetFieldValue(PositionAndTime, Point) ;
    TEvaluateRhsGivenB(y, Point, dydx);
}
#endif  // TMAGFIELDEQUATION_H 
