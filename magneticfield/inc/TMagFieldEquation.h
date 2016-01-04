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

// #include <vector>
#include "base/Vector3D.h"

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
     static constexpr double fCof   = Constants::c_light;  //   / fieldUnits::meter ;

     // Expected constant value:
     // static constexpr double fCof    = Constants::c_light * fieldUnits::second /
     //     ( 1.0e9 * fieldUnits::meter * fieldUnits::meter );

     TMagFieldEquation(T_Field* pF) : GUVEquationOfMotion(pF) { fPtrField = pF; }

     TMagFieldEquation(const TMagFieldEquation& );
     ~TMagFieldEquation()  {}  // Was virtual - but now no inheritance

     TMagFieldEquation<Field,Size>* Clone() const;
     TMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool& safe);
     TMagFieldEquation<Field,Size>* CloneOrSafeSelf(bool* safe=0);
        // If class is thread safe, return self ptr.  Else return clone
     
     REALLY_INLINE  // inline __attribute__((always_inline))
     void GetFieldValue(const double Point[4],
                              double Value[]) const
     {
        fPtrField->T_Field::GetFieldValue(Point, Value);
     }

     T_Field* GetField() { return fPtrField; } 

     inline // REALLY_INLINE
     void RightHandSide(const double y[], /*double charge,*/ double dydx[] ) const;

     inline // REALLY_INLINE
     void RightHandSide(const double y[], /*double charge,*/ double dydx[],
                        vecgeom::Vector3D<float> &BfieldVec) const;

     inline // REALLY_INLINE
     void RightHandSide(const             double  y[],
                        const vecgeom::Vector3D<double> &Position,   // Copy of y[0], y[1], y[2] - for optim
                                        /*double  charge,*/
                                          double  dydx[],
                        vecgeom::Vector3D<float> &BfieldVec) const;

     REALLY_INLINE
     void TEvaluateRhsGivenB( const double y[],
                              const vecgeom::Vector3D<float> B,  // Was double B[3],
                                //  double charge,
                             double dydx[] ) const;

     // virtual
     void EvaluateRhsGivenB( const double y[],
                             const vecgeom::Vector3D<float> B,  // Was const double B[3],
                              //   double charge,
                             double dydx[] ) const
     { TEvaluateRhsGivenB( y, B, /*charge,*/ dydx); }

     REALLY_INLINE
     void FieldFromY(const double y[], double Bfield[] ) const;

     REALLY_INLINE
     void FieldFromY(const double y[], vecgeom::Vector3D<float> &Bfield ) const;

     REALLY_INLINE
     void PrintInputFieldAndDyDx(const double y[], /* double charge, */ double dydx[] ) const;

     REALLY_INLINE
     void PrintAll(  double const  y[],
                     double const  B[3],
                     double        charge,
                     double        cof,                     
                     double const  dydx[]  ) const;
     
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
   TMagFieldEquation<Field,Size>
   ::TMagFieldEquation(const TMagFieldEquation& right)
   :
      GUVEquationOfMotion( (GUVField*) 0 ),
      fPtrField( right.fPtrField->CloneOrSafeSelf( (bool *)0 ) )
      // fPtrField( new Field(right.fPtrField) )
{
   // G4bool threadSafe;
   // fPtrField = right.fPtrField->CloneOrSafeSelf( &threadSafe );

   // std::cout <<  "TMagFieldEquation - copy constructor called." << std::endl;
   GUVEquationOfMotion::SetFieldObj( fPtrField ); //  Also stored in base class ... for now
}

template
<class Field, unsigned int Size>
   TMagFieldEquation<Field,Size>*
   TMagFieldEquation<Field,Size>
   ::Clone() const
{
   // bool safe= false;  // Field* pField= fPtrField->CloneOrSafeSelf(safe);
   Field* cloneField= fPtrField->Clone();   
   std::cerr << " #TMagFieldEquation<Field,Size>::Clone() called# " << std::endl;
   return new TMagFieldEquation( cloneField );
}

template
<class Field, unsigned int Size>
   TMagFieldEquation<Field,Size>*
   TMagFieldEquation<Field,Size>
   ::CloneOrSafeSelf(bool& safe)
{
   TMagFieldEquation<Field,Size>* equation;
   Field* pField=
      fPtrField->CloneOrSafeSelf(safe);
   // If Field does not have such a method:
   //  = new Field( fPtrField ); // Need copy constructor.
   //  safe= false;

   std::cerr << " #TMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool& safe) called# " << std::endl;

   // safe = safe && fClassSafe;
   // if( safe )  {  equation = this; }
   //    Can be introduced when Equation is thread safe -- no state
   //     --> For now the particle Charge is preventing this 23.11.2015
   // else {
      equation = new TMagFieldEquation( pField );
      safe= false;
   // }

   return equation;   
}

template
<class Field, unsigned int Size>
   TMagFieldEquation<Field,Size>*
   TMagFieldEquation<Field,Size>
   ::CloneOrSafeSelf(bool* pSafe)
{
   bool safeLocal;
   std::cerr << " #TMagFieldEquation<Field,Size>::CloneOrSafeSelf(bool* safe) called#" << std::endl;
   if( !pSafe ) pSafe= &safeLocal;
   auto equation= CloneOrSafeSelf( pSafe );
   return equation;   
}


template
<class Field, unsigned int Size>
REALLY_INLINE
   void  TMagFieldEquation<Field, Size>
   ::TEvaluateRhsGivenB( const double y[],
                         const vecgeom::Vector3D<float> Bfloat,  // Was const double B[3],
                           //  double charge,
                               double dydx[]  ) const
{
    const double charge = fParticleCharge; 
    // ThreeVectorD momentum( y[3], y[4], y[5]);
    double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    double inv_momentum_magnitude = 1. / std::sqrt( momentum_mag_square );
    // double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);
  
    double B[3]= { Bfloat[0], Bfloat[1], Bfloat[2] };

    double cof = charge*fCof*inv_momentum_magnitude;

//  printf("            B-field= %10.3f %10.3f %10.3f  mag= %10.3f %10.3f\n", B[0], B[1], B[2],
//        std::sqrt(Bmag2chk) , (double)Bfloat.Mag() );    
//  printf("           Momentum= %12.6g %12.6g %12.6g    mag= %12.7g \n", y[3], y[4], y[5], 1.0/inv_momentum_magnitude );
    
    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

    // PrintAll( y, B, charge, cof, dydx );
    return;
}


template
<class Field, unsigned int Size>
REALLY_INLINE
   void  TMagFieldEquation<Field, Size>
   ::PrintAll(  double const  y[],
                double const  B[3],
                double        charge,
                double        cof,
                double const  dydx[]  ) const
{
    using ThreeVectorD = vecgeom::Vector3D<double>;
   
    printf("Equation:  fCof = %8.4g  charge= %f cof= %10.5g   B-field= %f %f %f \n",
           fCof, charge, cof, B[0], B[1], B[2] );
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

template
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size>
   ::FieldFromY(const double y[],  double Bfield[3] ) const
{
   // double  PositionAndTime[4] = { y[0], y[1], y[2], 0.0 } 
   // PositionAndTime[3] = y[7];  //  --> extersion?
   // GetFieldValue(PositionAndTime, Bfield) ;
   GetFieldValue(y, Bfield) ;
}

template
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size>
   ::FieldFromY(const double y[],  vecgeom::Vector3D<float> &Bfield ) const
{
   // double  Bfield[3];  //G4maximum_number_of_field_components];
   // double  PositionAndTime[4] = { y[0], y[1], y[2], 0.0 } 
   // PositionAndTime[3] = y[7];  //  --> Time / extension?
   vecgeom::Vector3D<double> Position( y[0], y[1], y[2] );

   fPtrField->T_Field::GetFieldValue( Position, Bfield );
}


template
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size>
   ::RightHandSide(const double y[], /* double charge, */ double dydx[] ) const
{
    // double  BfieldArr[3];  //G4maximum_number_of_field_components];
    // FieldFromY( y, BfieldArr );
    // TEvaluateRhsGivenB( y, BfieldArr, /*charge,*/ dydx );

    vecgeom::Vector3D<float> BfieldVec;

    FieldFromY( y, BfieldVec );
    TEvaluateRhsGivenB( y, BfieldVec, /*charge,*/ dydx );
}

template
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size>
   ::RightHandSide(const double  y[],
                      /* double  charge, */
                         double  dydx[],
       vecgeom::Vector3D<float> &BfieldVec ) const
{
    FieldFromY( y, BfieldVec );
    TEvaluateRhsGivenB( y, BfieldVec, dydx );
}

template
<class Field, unsigned int Size>
REALLY_INLINE
void
TMagFieldEquation<Field,Size>
   ::RightHandSide(const                   double  y[],
                   const vecgeom::Vector3D<double> &Position,
                                         /*double  charge,*/
                                           double  dydx[],
                         vecgeom::Vector3D<float> &BfieldVec) const
{
   fPtrField->T_Field::GetFieldValue( Position, BfieldVec );
   TEvaluateRhsGivenB( y, BfieldVec, dydx );
}



#include <iostream>   // For debuging only

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
#endif  // TMAGFIELDEQUATION_H
