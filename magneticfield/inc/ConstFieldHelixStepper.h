/*
 * ConstFieldHelixStepper.h
 *
 *  Created on: January 29, 2016
 *      Author: J. Apostolakis
 *      
 *  Started from ConstBzFieldHelixStepper by S. Wenzel
 */

#ifndef CONSTVECFIELDHELIXSTEPPER_H_
#define CONSTVECFIELDHELIXSTEPPER_H_

#include <base/Vector3D.h>
#include <base/SOA3D.h>
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
  class ConstFieldHelixStepper
  {

    public:
      VECCORE_ATT_HOST_DEVICE
      ConstFieldHelixStepper( double Bx, double By, double Bz );

      VECCORE_ATT_HOST_DEVICE
      ConstFieldHelixStepper( double Bfield[3] );

      VECCORE_ATT_HOST_DEVICE
      ConstFieldHelixStepper( vecgeom::Vector3D<double> Bfield );      

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
       template<typename Vector3D, typename Real_v, typename BaseIType>
       GEANT_FORCE_INLINE
       VECCORE_ATT_HOST_DEVICE
          // inline __attribute__((always_inline))
       void DoStep( Real_v const & /*posx*/, Real_v const & /*posy*/, Real_v const & /*posz*/,
                    Real_v const & /*dirx*/, Real_v const & /*diry*/, Real_v const & /*dirz*/,
                    BaseIType const & /*charge*/, Real_v const & /*momentum*/, Real_v const & /*step*/,
                    Real_v & /*newsposx*/, Real_v  & /*newposy*/, Real_v  & /*newposz*/,
                    Real_v & /*newdirx*/, Real_v  & /*newdiry*/, Real_v  & /*newdirz*/
          ) const ;

       /**
        * basket version of dostep
        * version that takes plain arrays as input; suited for first-gen GeantV
        */
       template <typename Real_v>
       GEANT_FORCE_INLINE
       void
       DoStepArr( double const * posx, double const * posy, double const * posz,
                  double const * dirx, double const * diry, double const * dirz,
                  int    const * charge, double const * momentum, double const * step,
                  double * newsposx, double * newposy, double * newposz,
                  double * newdirx, double * newdiry, double * newdirz,
                  int np
          ) const ;

       /******
       template <typename Real_v>
       GEANT_FORCE_INLINE
       void
          DoStepSOA3D( const vecgeom::SOA3D<Real_v> & position,
                       const vecgeom::SOA3D<Real_v> & direction,
                       double const * momentumMag,
                       int    const * charge,
                       double const * step,
                             vecgeom::SOA3D<Real_v> & positionEnd,
                             vecgeom::SOA3D<Real_v> & directionEnd, 
                       int           numTracks
          ) const ;
        ******/
       
        // in future will offer versions that take containers as input

        /**
         * this function propagates the track along the helix solution by a step
         * input: current position, current direction, some particle properties
         * output: new position, new direction of particle
         */
       template<typename Vector3D_t, typename Real_v, typename BaseIType>
          GEANT_FORCE_INLINE
          void DoStep( Vector3D_t  const & position,
                       Vector3D_t  const & direction,
                       BaseIType const & charge,
                       Real_v    const & momentum,
                       Real_v    const & step,
                       Vector3D_t        & newPosition,
                       Vector3D_t        & newDirection
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
void ConstFieldHelixStepper::CalculateDerived()
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
ConstFieldHelixStepper::ConstFieldHelixStepper( double Bx, double By, double Bz )
   :  fBx(Bx), fBy(By), fBz(Bz)
  {
     CalculateDerived();
  }

inline
ConstFieldHelixStepper::ConstFieldHelixStepper( double B[3] )
   :  fBx(B[0]), fBy(B[1]), fBz(B[2])
  {
     CalculateDerived();
  }

inline
VECCORE_ATT_HOST_DEVICE   
ConstFieldHelixStepper::ConstFieldHelixStepper( vecgeom::Vector3D<double> Bfield )
   : fBx( Bfield.x() ), fBy( Bfield.y() ), fBz( Bfield.z() )
  {
     CalculateDerived();   
  }

  /**
   * this function propagates the track along the "helix-solution" by a step step
   * input: current position (x0, y0, z0), current direction ( dirX0, dirY0, dirZ0 ), some particle properties
   * output: new position, new direction of particle
   */
template<typename Vector3D_t, typename BaseDType, typename BaseIType>
   inline
   __attribute__((always_inline))
   void ConstFieldHelixStepper::DoStep(
               BaseDType const & x0, BaseDType const & y0, BaseDType const & z0,
               BaseDType const & dirX0, BaseDType const & dirY0, BaseDType const & dirZ0,
               BaseIType const & charge, BaseDType const & momentum, BaseDType const & step,
               BaseDType & x, BaseDType & y, BaseDType & z,
               BaseDType & dx, BaseDType & dy, BaseDType & dz
             ) const
  {
     Vector3D_t startPosition = { x0, y0, z0 };
     Vector3D_t startDirection( dirX0, dirY0, dirZ0 );
     Vector3D_t   endPosition,   endDirection;

     // startPosition.Set( x0, y0, z0);
     // startDirection.Set( dirX0, dirY0, dirZ0);

     DoStep( startPosition, startDirection, charge, momentum, step,
             endPosition, endDirection);
     x= endPosition.x();
     y= endPosition.y();
     z= endPosition.z();
     dx= endDirection.x();
     dy= endDirection.y();
     dz= endDirection.z();
  }

  template<typename Vector3D_t, typename BaseDType, typename BaseIType>
  inline
  __attribute__((always_inline))
  void ConstFieldHelixStepper::DoStep( Vector3D_t  const & startPosition,
                                       Vector3D_t  const & startDirection,
                                       BaseIType const & charge,
                                       BaseDType const & momentum,
                                       BaseDType const & step,
                                       Vector3D_t        & endPosition,
                                       Vector3D_t        & endDirection
     ) const
  {
      const double kB2C_local = -0.299792458e-3;
      const double kSmall = 1.E-30;
      using vecCore::math::Max;
      using vecCore::math::Sin;
      using vecCore::math::Cos;
      using vecCore::math::Abs;      
      // could do a fast square root here

      // BaseDType dt = sqrt((dx0*dx0) + (dy0*dy0)) + kSmall;

      // assert( std::abs( startDirection.Mag2() - 1.0 ) < 1.0e-6 );

      Vector3D_t  dir1Field( fUnitX, fUnitY, fUnitZ );
      BaseDType UVdotUB = startDirection.Dot(dir1Field);   //  Limit cases 0.0 and 1.0
      BaseDType dt2   = Max( startDirection.Mag2() - UVdotUB * UVdotUB, BaseDType(0.0) );
      BaseDType sinVB = sqrt( dt2 ) + kSmall;
 
      // BaseDType invnorm = 1. / sinVB;

      // radius has sign and determines the sense of rotation
      BaseDType R = momentum*sinVB/((kB2C_local*BaseDType(charge))*(fBmag));

      Vector3D_t  restVelX = startDirection - UVdotUB * dir1Field;

      // Vector3D_t  dirVelX( 0.0, 0.0, 0.0 );            // OK if it is zero - ie. dir // B
      // if( restVelX.Mag2() > 0.0 ) dirVelX = restVelX.Unit();
      Vector3D_t  dirVelX  = restVelX.Unit();             // Unit must cope with 0 length !!
      Vector3D_t  dirCrossVB = dirVelX.Cross(dir1Field);  // OK if it is zero 
      // Vector3D_t  dirCrossVB= restVelX.Cross(dir1Field);  // OK if it is zero      

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
      assert ( vecCore::MaskFull( Abs( dir1Field.Dot(  dirVelX   ) ) < 1.e-6 ) );
      assert ( vecCore::MaskFull( Abs( dirVelX.Dot(    dirCrossVB) ) < 1.e-6 ) );
      assert ( vecCore::MaskFull( Abs( dirCrossVB.Dot( dir1Field ) ) < 1.e-6 ) );
      
      BaseDType phi = - step * BaseDType(charge) * fBmag * kB2C_local / momentum;

      BaseDType cosphi = Cos(phi);
      BaseDType sinphi = Sin(phi);
      // sincos(phi, &sinphi, &cosphi);
      // SinCos(phi, &sinphi, &cosphi);  // Opportunity for new 'efficient' method !? 
      
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

//  SW: commented out due to explicit Vc dependence and since it is not currently used
//       leaving the code here to show how one would dispatch to the kernel with Vc
#define _R_ __restrict__

  template <typename Real_v>
  void
  ConstFieldHelixStepper::DoStepArr(
                        double const * _R_ posx, double const * _R_ posy, double const * _R_ posz,
                        double const * _R_ dirx, double const * _R_ diry, double const * _R_ dirz,
                        int const * _R_ charge, double const * _R_ momentum, double const * _R_ step,
                        double * _R_ newposx, double * _R_ newposy, double * _R_ newposz,
                        double * _R_ newdirx, double * _R_ newdiry, double * _R_ newdirz,
                        int np
                     ) const
   {
       const int vectorSize= vecgeom::kVectorSize;
       using vecCore::Load;
       using vecCore::Store;
       using vecCore::Set;

       int i;
       for ( i=0; i < np ; i+= vectorSize )
       {
            // results cannot not be temporaries
            //    Vector3D<Real_v> newPosition, newDirection;
            Real_v oldPosx_v, oldPosy_v, oldPosz_v,
                   oldDirx_v, oldDiry_v, oldDirz_v;
            Real_v newposx_v, newposy_v, newposz_v,
                   newdirx_v, newdiry_v, newdirz_v;
            Real_v momentum_v, stepSz_v;
            
            vecCore::Index<Real_v> chargeVint;
            Load( chargeVint, &charge[i] );  // ==> Does not work for index Type !!
            for( size_t j = 0; j < vectorSize; ++j )
               Set( chargeVint, j, charge[i+j] );
               
            Load( oldPosx_v, &posx[i] );
            Load( oldPosy_v,  &posy[i] );
            Load( oldPosz_v,  &posz[i] );
            Load( oldDirx_v,  &dirx[i] );
            Load( oldDiry_v,  &diry[i] );
            Load( oldDirz_v,  &dirz[i] );
            Load( momentum_v, &momentum[i] );
            Load( stepSz_v,   &step[i] );
            
            DoStep<vecgeom::Vector3D<Real_v>, Real_v, vecCore::Index<Real_v>>
               ( oldPosx_v, oldPosy_v, oldPosz_v,
                    // Real_v(posx[i]), Real_v(posy[i]), Real_v(posz[i]),
                 oldDirx_v, oldDiry_v, oldDirz_v,
                    // Real_v(dirx[i]), Real_v(diry[i]), Real_v(dirz[i]),
                 chargeVint, // Geant::Int_v(charge[i]),
                 momentum_v, // Real_v(momentum[i]),
                 stepSz_v,   // Real_v(step[i]),
                 newposx_v, newposy_v, newposz_v,
                 newdirx_v, newdiry_v, newdirz_v
               );
            // write results
            Store(newposx_v, &newposx[i]);
            Store(newposy_v, &newposy[i]);
            Store(newposz_v, &newposz[i]);
            Store(newdirx_v, &newdirx[i]);
            Store(newdiry_v, &newdiry[i]);
            Store(newdirz_v, &newdirz[i]);
       }

       // tail part
       for (; i < np ; i++ )
          DoStep<vecgeom::Vector3D<double>, double, vecCore::Index<double>>
             ( posx[i], posy[i], posz[i],
                  dirx[i], diry[i], dirz[i],
                  charge[i],
                  momentum[i],
                  step[i],
                  newposx[i],
                  newposy[i],
                  newposz[i],
                  newdirx[i],
                  newdiry[i],
                  newdirz[i]
             );
   }

 /********  Attempt to create method directly for SOA3D 
  template <typename Real_v>
     GEANT_FORCE_INLINE
     void ConstFieldHelixStepper::
     DoStepSOA3D( const vecgeom::SOA3D<Real_v> & position,
                  const vecgeom::SOA3D<Real_v> & direction,
                  double const * momentumMag,
                  int    const * charge,
                  double const * step,
                  vecgeom::SOA3D<Real_v> & positionEnd,
                  vecgeom::SOA3D<Real_v> & directionEnd, 
                       int           numTracks
        ) const
  {
     // template <typename Real_vecTp> using Vector3D = vecgeom::Vector3D<Real_vecTp>;
     
     // Use the values in the SOA3D directly - without minimum of copying
     const int vectorSize= vecgeom::VectorSize;
     int i;
     for ( i=0; i < numTracks ; i+= vectorSize )
     {
        // results cannot not be temporaries
        vecgeom::
        Vector3D<Real_v> newPosition, newDirection;           
        
        vecCore::Index<Real_v> chargeVint;
        Load( chargeVint, charge[i] );

        vecgeom::        
        Vector3D<Real_v> positionVec= position[i]; // ( position.fX[i],  position.fY[i],  position.fZ[i] );
        // Vector3D<Real_v> directionVec( direction.fX[i], direction.fY[i], direction.fZ[i] );
        
            // Cannot masquarade for the end values !!!
            // Vector3D<Real_v> endPositionVec(   position.fX+i,  position.fY+i,  position.fZ+i );
            // Vector3D<Real_v> endDirectionVec( direction.fX+i, direction.fY+i, direction.fZ+i );
            Real_v momentumV = Real_v( ); 
            DoStep( position[i], // positionVec,
                    direction[i], // directionVec,
                    chargeVint,
                    Real_v( &momentum[i] ),
                    Real_v( &step[i] ),
                    endPositionVec,
                    endDirectionVec
                   );
            positionEnd[i] = endPositionVec;                                         
            directionEnd[i] = endDirectionVec;                                           
            // StoreInSOA3D(  positionEnd, i,  endPositionVec );
            // StoreInSOA3D( directionEnd, i, endDirectionVec );
       }

       // tail part
       for (; i < numTracks ; i++ )
          DoStep<double>( posx[i], posy[i], posz[i],
                  dirx[i], diry[i], dirz[i],
                  charge[i],
                  momentum[i],
                  step[i],
                  newposx[i],
                  newposy[i],
                  newposz[i],
                  newdirx[i],
                  newdiry[i],
                  newdirz[i]
             );     
  }
  *************/
  
} // inline namespace GEANT_IMPL_NAMESPACE
} // end geant namespace

#endif /* CONSTVECFIELDHELIXSTEPPER_H_ */
