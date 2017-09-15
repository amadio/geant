/*
 * ConstVecFieldHelixStepper.h
 *
 *  Created on: January 29, 2016
 *      Author: J. Apostolakis
 *      
 *  Started from ConstBzFieldHelixStepper by S. Wenzel
 */

#ifndef CONSTVECFIELDHELIXSTEPPER_H_
#define CONSTVECFIELDHELIXSTEPPER_H_

#include <base/Vector3D.h>
#include <Geant/Config.h>
#include <Geant/VectorTypes.h>

#include "GFldAuxFunctions.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
  /**
  * A very simple stepper treating the propagation of particles in a constant Bz magnetic field
  * ( neglecting energy loss of particle )
  * This class is roughly equivalent to TGeoHelix in ROOT
  */
  class ConstVecFieldHelixStepper
  {

    public:
      VECCORE_ATT_HOST_DEVICE
      ConstVecFieldHelixStepper( double Bx, double By, double Bz );

      VECCORE_ATT_HOST_DEVICE
      ConstVecFieldHelixStepper( double Bfield[3] );

      void SetBx( double Bx ){ fBx = Bx; CalculateDerived(); }
      void SetBy( double By ){ fBy = By; CalculateDerived(); }
      void SetBz( double Bz ){ fBz = Bz; CalculateDerived(); }
      double GetBz() const { return fBz; }

      /*
      template<typename RT, typename Vector3D>
      RT GetCurvature( Vector3D const & dir,
              double const charge, double const momentum ) const
      {
        if (charge == 0) return RT(0.);
        return abs( kB2C * fBz * dir.FastInverseScaledXYLength( momentum ) );
      }
      */

      /**
       * this function propagates the track along the helix solution by a step
       * input: current position, current direction, some particle properties
       * output: new position, new direction of particle
       */
       template<typename Vector3D, typename BaseType, typename BaseIType>
       inline
       __attribute__((always_inline))
       VECCORE_ATT_HOST_DEVICE
       void DoStep( BaseType const & /*posx*/, BaseType const & /*posy*/, BaseType const & /*posz*/,
                    BaseType const & /*dirx*/, BaseType const & /*diry*/, BaseType const & /*dirz*/,
                    BaseIType const & /*charge*/, BaseType const & /*momentum*/, BaseType const & /*step*/,
                    BaseType & /*newsposx*/, BaseType  & /*newposy*/, BaseType  & /*newposz*/,
                    BaseType & /*newdirx*/, BaseType  & /*newdiry*/, BaseType  & /*newdirz*/
          ) const ;

       /**
        * basket version of dostep
        * version that takes plain arrays as input; suited for current Geant-V
        *
        * SW: for the moment (12.5.2015) commenting this out as not used 
        *
        */
       //void DoStep_v( double const * /*posx*/, double const * /*posy*/, double const * /*posz*/,
       //               double const * /*dirx*/, double const * /*diry*/, double const * /*dirz*/,
       //              int const * /*charge*/, double const * /*momentum*/, double const * /*step*/,
       //              double * /*newsposx*/, double * /*newposy*/, double * /*newposz*/,
       //              double * /*newdirx*/, double * /*newdiry*/, double * /*newdirz*/,
       //              int np
       //           ) const ;

        // in future will offer versions that take containers as input

        /**
         * this function propagates the track along the helix solution by a step
         * input: current position, current direction, some particle properties
         * output: new position, new direction of particle
         */
         template<typename Vector3D, typename BaseType, typename BaseIType>
         inline
         __attribute__((always_inline))
         void DoStep( Vector3D const & pos, Vector3D const & dir,
                           BaseIType const & charge, BaseType const &  momentum,
                           BaseType const & step,
                           Vector3D & newpos,
                           Vector3D & newdir
            ) const;

    protected:
      void CalculateDerived();

    private:
      double fBx, fBy, fBz;
      // Auxilary members - calculated from above - cached for speed, code simplicity
      double fBmag;
      double fUnitX, fUnitY, fUnitZ;
  }; // end class declaration

inline 
void ConstVecFieldHelixStepper::CalculateDerived()
  {
     fBmag = std::sqrt( fBx * fBx + fBy * fBy + fBz * fBz );
     if ( fBmag > 0.0 ) {
        double inv_mag=  1.0 / fBmag;
        fUnitX = inv_mag * fBx;
        fUnitY = inv_mag * fBy;
        fUnitZ = inv_mag * fBz;
     } else {
        // Arbitrary unit vector
        fUnitX = 0.0;
        fUnitY = 0.0; 
        fUnitZ = 1.0;
     }
  }

inline
ConstVecFieldHelixStepper::ConstVecFieldHelixStepper( double Bx, double By, double Bz )
   :  fBx(Bx), fBy(By), fBz(Bz)
  {
     CalculateDerived();
  }

inline
ConstVecFieldHelixStepper::ConstVecFieldHelixStepper( double B[3] )
   :  fBx(B[0]), fBy(B[1]), fBz(B[2])
  {
     CalculateDerived();
  }

  /**
   * this function propagates the track along the "helix-solution" by a step step
   * input: current position (x0, y0, z0), current direction ( dx0, dy0, dz0 ), some particle properties
   * output: new position, new direction of particle
   */
template<typename Vector3D, typename BaseDType, typename BaseIType>
   inline
   __attribute__((always_inline))
   void ConstVecFieldHelixStepper::DoStep(
               BaseDType const & x0, BaseDType const & y0, BaseDType const & z0,
               BaseDType const & dx0, BaseDType const & dy0, BaseDType const & dz0,
               BaseIType const & charge, BaseDType const & momentum, BaseDType const & step,
               BaseDType & x, BaseDType & y, BaseDType & z,
               BaseDType & dx, BaseDType & dy, BaseDType & dz
             ) const
  {
     Vector3D startPosition( x0, y0, z0);
     Vector3D startDirection( dx0, dy0, dz0);
     Vector3D endPosition, endDirection;

     DoStep( startPosition, startDirection, charge, momentum, step,
             endPosition, endDirection);
     x= endPosition.x();
     y= endPosition.y();
     z= endPosition.z();
     dx= endDirection.x();
     dy= endDirection.y();
     dz= endDirection.z();
  }

  template<typename Vector3D, typename BaseDType, typename BaseIType>
  inline
  __attribute__((always_inline))
  void ConstVecFieldHelixStepper::DoStep( Vector3D  const & startPosition,
                                          Vector3D  const & startDirection,
                                          BaseIType const & charge,
                                          BaseDType const & momentum,
                                          BaseDType const & step,
                                          Vector3D        & endPosition,
                                          Vector3D        & endDirection
     ) const
  {
      const double kB2C_local = -0.299792458e-3;
      const double kSmall = 1.E-30;
      // could do a fast square root here

      // BaseDType dt = sqrt((dx0*dx0) + (dy0*dy0)) + kSmall;

      assert( std::abs( startDirection.Mag2() - 1.0 ) < 1.0e-6 );

      Vector3D  dir1Field( fUnitX, fUnitY, fUnitZ );
      BaseDType UVdotUB = startDirection.Dot(dir1Field);   //  Limit cases 0.0 and 1.0
      BaseDType dt2   = std::max( startDirection.Mag2() - UVdotUB * UVdotUB, 0.0 );
      BaseDType sinVB = sqrt( dt2 ) + kSmall;
 
      // BaseDType invnorm = 1. / sinVB;

      // radius has sign and determines the sense of rotation
      BaseDType R = momentum*sinVB/((kB2C_local*BaseDType(charge))*(fBmag));

      Vector3D  restVelX = startDirection - UVdotUB * dir1Field;

      // Vector3D  dirVelX( 0.0, 0.0, 0.0 );            // OK if it is zero - ie. dir // B
      // if( restVelX.Mag2() > 0.0 ) dirVelX = restVelX.Unit();
      Vector3D  dirVelX  = restVelX.Unit();             // Unit must cope with 0 length !!
      Vector3D  dirCrossVB = dirVelX.Cross(dir1Field);  // OK if it is zero 
      // Vector3D  dirCrossVB= restVelX.Cross(dir1Field);  // OK if it is zero      

      /***
      printf("\n");
      printf("CVFHS> dir-1  B-fld  = %f %f %f   mag-1= %g \n", dir1Field.x(), dir1Field.y(), dir1Field.z(),
             dir1Field.Mag()-1.0 );
      printf("CVFHS> dir-2  VelX   = %f %f %f   mag-1= %g \n", dirVelX.x(), dirVelX.y(), dirVelX.z(),
             dirVelX.Mag()-1.0 );
      printf("CVFHS> dir-3: CrossVB= %f %f %f   mag-1= %g \n", dirCrossVB.x(), dirCrossVB.y(), dirCrossVB.z(),
             dirCrossVB.Mag()-1.0 );
      // dirCrossVB = dirCrossVB.Unit();
      printf("CVFHS> Dot products   d1.d2= %g   d2.d3= %g  d3.d1= %g \n",
             dir1Field.Dot(dirVelX), dirVelX.Dot( dirCrossVB), dirCrossVB.Dot(dir1Field) );
       ***/
      assert ( std::fabs( dir1Field.Dot(  dirVelX   ) ) < 1.e-6 );
      assert ( std::fabs( dirVelX.Dot(    dirCrossVB) ) < 1.e-6 );      
      assert ( std::fabs( dirCrossVB.Dot( dir1Field ) ) < 1.e-6 );
      
      BaseDType phi = - step * BaseDType(charge) * fBmag * kB2C_local / momentum;

      BaseDType cosphi;
      BaseDType sinphi;
      sincos(phi, &sinphi, &cosphi);

      endPosition = startPosition + R * ( cosphi - 1 ) * dirCrossVB
                    - R * sinphi * dirVelX
                    + step * UVdotUB * dir1Field;   //   'Drift' along field direction

      // dx = dx0 * cosphi - sinphi * dy0;
      // dy = dy0 * cosphi + sinphi * dx0;
      // dz = dz0;
      // printf(" phi= %f, sin(phi)= %f , sin(V,B)= %f\n", phi, sinphi, sinVB );
      endDirection = UVdotUB * dir1Field + cosphi * sinVB * dirVelX + sinphi * sinVB * dirCrossVB;
  }
 
  /**
   * basket version of dostep
   */
  /*
   SW: commented out due to explicit Vc dependence and since it is not currently used
       leaving the code here to show how one would dispatch to the kernel with Vc
#define _R_ __restrict__
  void ConstVecFieldHelixStepper::DoStep_v(
                        double const * _R_ posx, double const * _R_ posy, double const * _R_ posz,
                        double const * _R_ dirx, double const * _R_ diry, double const * _R_ dirz,
                        int const * _R_ charge, double const * _R_ momentum, double const * _R_ step,
                        double * _R_ newposx, double * _R_ newposy, double * _R_ newposz,
                        double * _R_ newdirx, double * _R_ newdiry, double * _R_ newdirz,
                        int np
                     ) const
   {
       // we have choice here: ( try autovectorization: )

//#pragma ivdep
//      for (int i=0;i<np;++i){
//            DoStep( posx[i], posy[i], posz[i], dirx[i], diry[i], dirz[i],
//                    charge[i], momentum[i], step[i],
//                    newposx[i], newposy[i], newposz[i],
//                    newdirx[i], newdiry[i], newdirz[i]
//                  );
//       }

       // alternative loop with Vc:

       for (int i=0;i<np;i+=Geant::kVecLenD)
       {
            // results cannot not be temporaries
            Double_v newposx_v, newposy_v, newposz_v,
                       newdirx_v, newdiry_v,newdirz_v;
            DoStep( Double_v(posx[i]), Double_v(posy[i]), Double_v(posz[i]),
                    Double_v(dirx[i]), Double_v(diry[i]), Double_v(dirz[i]),
                    Geant::Int_v(charge[i]), Double_v(momentum[i]), Double_v(step[i]),
                    newposx_v,
                    newposy_v,
                    newposz_v,
                    newdirx_v,
                    newdiry_v,
                    newdirz_v
                   );
            // write results
            newposx_v.store(&newposx[i]);
            newposy_v.store(&newposy[i]);
            newposz_v.store(&newposz[i]);
            newdirx_v.store(&newdirx[i]);
            newdiry_v.store(&newdiry[i]);
            newdirz_v.store(&newdirz[i]);
       }
       // tail part: tobedone
   }
  */

   //TODO: above stepper is tailored/specialized to B=(0,0,Bz) in the global frame of reference
   // might need to provide more general class in which the constant field has arbitrary direction

} // inline namespace GEANT_IMPL_NAMESPACE
} // end geant namespace

#endif /* CONSTVECFIELDHELIXSTEPPER_H_ */
