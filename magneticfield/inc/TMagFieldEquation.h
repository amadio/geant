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

#include "Units.h"
#include "Constants.h"
//  Update to GeantV units ASAP

template 
<class Field, unsigned int Size>
class TMagFieldEquation : public GUVEquationOfMotion
{
   public:
     typedef Field T_Field;
     static const unsigned int  N   = Size;
     static constexpr double fCof    = fieldUnits::eplus * Constants::c_light ;
   
     TMagFieldEquation(T_Field* pF) : GUVEquationOfMotion(pF) { fPtrField = pF; }
     ~TMagFieldEquation()  {}  // Was virtual - but now no inheritance

     REALLY_INLINE  // inline __attribute__((always_inline))     
     void GetFieldValue(const double Point[4],
                              double Value[]) const
     {
        fPtrField->T_Field::GetFieldValue(Point, Value);
     }

     inline // REALLY_INLINE
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

     REALLY_INLINE
     void FieldFromY(const double y[], /* double charge, */ double Bfield[] ) const;
   
     REALLY_INLINE
     void PrintInputFieldAndDyDx(const double y[], /* double charge, */ double dydx[] ) const;

     REALLY_INLINE
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
   ::FieldFromY(const double y[], /* double charge, */ double Bfield[3] ) const
{
    // double  Bfield[3];  //G4maximum_number_of_field_components];
    double  PositionAndTime[4];
    PositionAndTime[0] = y[0];
    PositionAndTime[1] = y[1];
    PositionAndTime[2] = y[2];
    PositionAndTime[3] = 0;        // Time
    // PositionAndTime[3] = y[7];  //  --> extersion?
    
    GetFieldValue(PositionAndTime, Bfield) ;
}

template 
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size> 
   ::RightHandSide(const double y[], /* double charge, */ double dydx[] ) const
{
    double  Bfield[3];  //G4maximum_number_of_field_components];

#if 0
    double  PositionAndTime[4];
    PositionAndTime[0] = y[0];
    PositionAndTime[1] = y[1];
    PositionAndTime[2] = y[2];
    PositionAndTime[3] = 0;        // Time
    // PositionAndTime[3] = y[7];  //  --> extersion?
    GetFieldValue(PositionAndTime, Bfield) ;
#endif
    FieldFromY( y, Bfield );
    
    TEvaluateRhsGivenB(y, Bfield, dydx);
}


#include <iostream>   // For debuging only
using std::cout;
using std::endl;

template 
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size> 
   ::PrintInputFieldAndDyDx(const double y[], /* double charge, */ double dydx[] ) const
{

    RightHandSide(y, dydx);

    // Obtain the field value 
    double  Bfield[3];  //G4maximum_number_of_field_components];
    FieldFromY( y, Bfield );
    TEvaluateRhsGivenB(y, Bfield, dydx);
    
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
