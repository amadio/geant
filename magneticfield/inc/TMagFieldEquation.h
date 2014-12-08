// Approach is derived from the Geant4 class G4MagFieldEquation
// 

// #include "G4ChargeState.hh"
// #include "G4Mag_UsualEqRhs.hh"
// #include "alt_sqrt.h"

//  Ensure that equation Right Hand Side is inlined - may be compiler dependend
#define INLINERHS 1

template 
<class Field, size_t Size>
class TMagFieldEquation // : public G4Mag_UsualEqRhs
{
    public:

        typedef Field T_Field;
        static const size_t N = Size;

        TMagFieldEquation(T_Field* f) { itsField = f; }
        ~TMagFieldEquation()  {}  // Was virtual - but now no inheritance

        inline __attribute__((always_inline)) 
        void GetFieldValue(const double Point[4],
                                 double Value[]) const
        {
            itsField->T_Field::GetFieldValue(Point, Value);
        }

   
  #ifdef INLINERHS
   // #pragma message "INLINING RHS"
        inline __attribute__((always_inline)) 
   //   #else
   // #pragma message "NOT INLINING RHS"
  #endif
        void RightHandSide(const double y[], double dydx[] ) const;


  //#ifdef INLINERHS
        inline __attribute__((always_inline))
  // #endif
        void TEvaluateRhsGivenB( const double y[],
                const double B[3],
                      double dydx[] ) const
        {
            double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
	    double inv_momentum_magnitude = 1./std::sqrt( momentum_mag_square);
	    //            double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 2);
            double cof = FCof()*inv_momentum_magnitude;

            dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
            dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
            dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

            dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
            dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
            dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

            return ;
        }

    private:
        enum { G4maximum_number_of_field_components = 24 };
        T_Field *itsField;
};


template 
<class Field, size_t Size>
void
#ifdef INLINERHS
       inline __attribute__((always_inline)) 
#else
#endif
TMagFieldEquation<Field,Size>::RightHandSide(const double y[], double dydx[] ) const
{
	  double Point[4];  //G4maximum_number_of_field_components]; 
	  double  PositionAndTime[3];
	  PositionAndTime[0] = y[0];
	  PositionAndTime[1] = y[1];
	  PositionAndTime[2] = y[2];
	  // PositionAndTime[3] = y[7];    // Tim
	  GetFieldValue(PositionAndTime, Point) ;
	  TEvaluateRhsGivenB(y, Point, dydx);
}
