//
// class TemplateFieldTrack
//
// Class description:
//
// Data structure bringing together a magnetic track's state.
// (position, momentum direction & modulus, energy, spin, ... )
// Uses/abilities:
//  - does not maintain any relationship between its data (eg E/p)
//  - for use in Runge-Kutta solver (now in passing it the values).
//
// History:
//
// Adaptations of template interface: J. Apostolakis,     Nov 2017
// First templated version:  Ananya, Feb/March 2016
//     ( commit 95e1316bcc156a04c876d6ea0fc9e60a15eeac4f )
// 
// Adapted from G4MagInt_Drv class of Geant4 (G4MagIntegratorDriver)
//
// History  - derived from G4FieldTrack
// - First version: Dec 9, 2014 John Apostolakis
// - First templated version:  Ananya Jan/Feb 2017, (TemplateGUFieldTrack)
//
// - Contributors: Ananya, J.Apostolakis                    2015-2017
// -----------------------------------------------------------------------

#ifndef TemplateFieldTrack_HH
#define TemplateFieldTrack_HH

#include "base/Vector3D.h"   // VecGeom/base/Vector3D.h 

template <class Real_v>
class  TemplateFieldTrack
{

   public:  // with description

     using ThreeVectorSIMD= vecgeom::Vector3D<Real_v>;

     TemplateFieldTrack(const ThreeVectorSIMD& pPosition, 
                        const ThreeVectorSIMD& pMomentum,
                         // double             restMass_c2,
                         // Real_v             charge,
                        Real_v  laboratoryTimeOfFlight= 0.0,
                        Real_v            curve_length= 0.0); 

     TemplateFieldTrack( const TemplateFieldTrack<Real_v>&   pFieldTrack ); 
     TemplateFieldTrack( char );   //  Almost default constructor

     ~TemplateFieldTrack();
       // End of preferred Constructors / Destructor 

     inline void
     UpdateState( const ThreeVectorSIMD& pPosition, 
                        Real_v     LaboratoryTimeOfFlight,
                  const ThreeVectorSIMD& pMomentumDirection,
                        Real_v     momentum              ); 
        //  Update four-vectors for space/time and momentum/energy
        //    Also resets curve length.

     // void SetCharge(double charge) { fCharge= charge; }

     inline TemplateFieldTrack& operator = ( const TemplateFieldTrack<Real_v> & rStVec );
       // Assignment operator

     inline ThreeVectorSIMD  GetMomentum() const;
     inline ThreeVectorSIMD  GetPosition() const;
     inline ThreeVectorSIMD  GetMomentumDirection() const;
     inline Real_v     GetMomentumMag() const;
     inline Real_v     GetCurveLength() const;
       // Distance along curve of point.

     // inline ThreeVectorSIMD  GetPolarization()   const; 
     // inline void         SetPolarization( const ThreeVectorSIMD& vecPol );

     inline Real_v     GetLabTimeOfFlight() const;
     inline Real_v     GetProperTimeOfFlight() const;
     // inline double       GetKineticEnergy() const;
       // Accessors.

     inline void SetPosition(ThreeVectorSIMD nPos); 
     inline void SetMomentum(ThreeVectorSIMD nMom);
       // Does change mom-dir too.

     inline void SetCurvePnt(const ThreeVectorSIMD& pPosition, 
                             const ThreeVectorSIMD& pMomentum,  
                                   Real_v      s_curve );

     // inline void SetMomentumDir(ThreeVectorSIMD nMomDir);
     // Does NOT change Momentum or Velocity Vector.

     // inline void SetRestMass(double Mass_c2) { fRestMass_c2= Mass_c2; }

       // Access
     // inline double GetCharge() const { return fCharge; } 
   
     inline void SetCurveLength(Real_v nCurve_s);
       // Distance along curve.
     inline void SetKineticEnergy(Real_v nEnergy);
       // Does not modify momentum.

     inline void SetLabTimeOfFlight(Real_v tofLab); 
     inline void SetProperTimeOfFlight(Real_v tofProper);
       //  Modifiers

   public: // without description

     static constexpr int ncompSVEC = 12;
       // Needed and should be used only for RK integration driver

     inline void DumpToArray(Real_v valArr[ncompSVEC]) const; 
            void LoadFromArray(const Real_v valArr[ncompSVEC],
                                     int      noVarsIntegrated );
     template <class Backend_>
     friend  std::ostream&
             operator<<( std::ostream& os, const TemplateFieldTrack<Backend_>& SixVec);

  private: // public: by A. ?
     Real_v  fPositionMomentum[6];
     
  private:
     Real_v  fDistanceAlongCurve;  // distance along curve of point
     Real_v  fMomentumMag;
     // Real_v  fKineticEnergy;
     // Real_v  fRestMass_c2;
     Real_v  fLabTimeOfFlight;
     Real_v  fProperTimeOfFlight;
     // ThreeVectorSIMD fPolarization;
     // ThreeVectorSIMD fMomentumDir;
     // Real_v  fInitialMomentumMag;  // At 'track' creation.
     // Real_v  fLastMomentumMag;     // From last Update (for checking.)

     // Real_v fCharge;
}; 

// #include "TemplateFieldTrack.icc"

//
// $Id: TemplateFieldTrack.icc 81175 2014-05-22 07:39:10Z gcosmo $
//


template <class Real_v>
inline
TemplateFieldTrack<Real_v>::
  TemplateFieldTrack( const TemplateFieldTrack<Real_v>&  rStVec  )
    : fDistanceAlongCurve( rStVec.fDistanceAlongCurve),
      fMomentumMag( rStVec.fMomentumMag ),
      // fKineticEnergy( rStVec.fKineticEnergy ),
      // fRestMass_c2( rStVec.fRestMass_c2),
      fLabTimeOfFlight( rStVec.fLabTimeOfFlight ), 
      fProperTimeOfFlight( rStVec.fProperTimeOfFlight ) //, 
      // fMomentumModulus( rStVec.fMomentumModulus ),
      // fPolarization( rStVec.fPolarization ), 
      // fMomentumDir( rStVec.fMomentumDir ), 
      // fCharge( rStVec.fCharge )
{

  //try auto-vectorization
  for (int i = 0; i < 6; ++i)
  {
    fPositionMomentum[i] = rStVec.fPositionMomentum[i];
  }
}

template <class Real_v>
inline
TemplateFieldTrack<Real_v>::~TemplateFieldTrack()
{
}

template <class Real_v>
inline void
TemplateFieldTrack<Real_v>::
  SetCurvePnt(const vecgeom::Vector3D<Real_v> &pPosition, 
              const vecgeom::Vector3D<Real_v> &pMomentum,  
                                      Real_v     s_curve )
{
  //try auto-vectorization
  for (int i = 0; i < 3; ++i)
  {
    fPositionMomentum[i]   = pPosition[i];
    fPositionMomentum[i+3] = pMomentum[i];
  }

  fMomentumMag= pMomentum.Mag();
  //Commented block below because seems to do nothing. If required, use a MaskedAssign : Ananya
/*  if( fMomentumMag > 0.0 )
  {
     // fMomentumDir = (1.0/fMomentumMag) * pMomentum;
  }*/
  fDistanceAlongCurve= s_curve;
} 


template <class Real_v>
inline
vecgeom::Vector3D<Real_v> 
TemplateFieldTrack<Real_v>::
  GetPosition() const
{
   vecgeom::Vector3D<Real_v> myPosition( fPositionMomentum[0], fPositionMomentum[1], fPositionMomentum[2] );
   return myPosition;
} 

template <class Real_v>
inline
vecgeom::Vector3D<Real_v> 
TemplateFieldTrack<Real_v>::
  GetMomentumDirection() const 
{
   typedef vecgeom::Vector3D<Real_v> ThreeVectorSIMD;
   Real_v inv_mag= 1.0 / fMomentumMag;
   return inv_mag * ThreeVectorSIMD( fPositionMomentum[3], fPositionMomentum[4], fPositionMomentum[5] );
} 

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  SetPosition( vecgeom::Vector3D<Real_v> pPosition) 
{
  //try auto-vectorization
  for (int i = 0; i < 3; ++i)
  {
    fPositionMomentum[i] = pPosition[i];
  }
} 

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  SetMomentum( vecgeom::Vector3D<Real_v> vMomentum) 
{
  fMomentumMag= vMomentum.Mag();

  // try auto-vectorization
  for (int i = 0; i < 3; ++i) {
    fPositionMomentum[i+3] = vMomentum[i];
  }
  // ELSE: 
  // fPositionMomentum[3] = vMomentum[0];
  // fPositionMomentum[4] = vMomentum[1];
  // fPositionMomentum[5] = vMomentum[2];
} 

template <class Real_v>
inline
Real_v 
TemplateFieldTrack<Real_v>::
  GetMomentumMag() const 
{
   return fMomentumMag;
} 

template <class Real_v>
inline
Real_v  
TemplateFieldTrack<Real_v>::
  GetCurveLength() const 
{
   return  fDistanceAlongCurve;  
}

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  SetCurveLength(Real_v nCurve_s)
{
   fDistanceAlongCurve= nCurve_s;
}

// inline Real_v TemplateFieldTrack<Real_v>::GetKineticEnergy() const
// { return fKineticEnergy; }

// inline void TemplateFieldTrack<Real_v>::SetKineticEnergy(Real_v newKinEnergy)
// {  fKineticEnergy=newKinEnergy; }

// inline ThreeVectorSIMD TemplateFieldTrack<Real_v>::GetPolarization() const
// { return fPolarization; }

// inline void TemplateFieldTrack<Real_v>::SetPolarization(const ThreeVectorSIMD& vecPlz)
// { fPolarization= vecPlz; }

template <class Real_v>
inline
Real_v 
TemplateFieldTrack<Real_v>::
  GetLabTimeOfFlight() const
{
   return fLabTimeOfFlight;
}

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  SetLabTimeOfFlight(Real_v nTOF)
{
   fLabTimeOfFlight=nTOF;
}

template <class Real_v>
inline
Real_v  
TemplateFieldTrack<Real_v>:: 
  GetProperTimeOfFlight() const
{
   return fProperTimeOfFlight;
}

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  SetProperTimeOfFlight(Real_v nTOF)
{
   fProperTimeOfFlight=nTOF;
}

template <class Real_v>
inline
vecgeom::Vector3D<Real_v> 
TemplateFieldTrack<Real_v>::
  GetMomentum() const 
{
   return ThreeVectorSIMD( fPositionMomentum[3], fPositionMomentum[4], fPositionMomentum[5] );
} 

// Dump values to array
//  
//   note that momentum direction is not saved 

template <class Real_v>
inline
void TemplateFieldTrack<Real_v>::
  DumpToArray(Real_v valArr[ncompSVEC] ) const
{

  //try auto-vectorization
  for (int i = 0; i < 6; ++i)
  {
    valArr[i] = fPositionMomentum[i];
  }


  ThreeVectorSIMD Momentum(valArr[3],valArr[4],valArr[5]);

  // double mass_in_Kg;
  // mass_in_Kg = fEnergy / velocity_mag_sq * (1-velocity_mag_sq/c_squared);
  // valArr[6]= mass_in_Kg;

  // The following components may or may not be integrated.
  // valArr[6]= fKineticEnergy; 

  // valArr[6]=fEnergy;  // When it is integrated over, do this ...
  valArr[7] = fLabTimeOfFlight;
  valArr[8] = fProperTimeOfFlight;
  // valArr[9]=fPolarization.x();
  // valArr[10]=fPolarization.y();
  // valArr[11]=fPolarization.z();
  // valArr[]=fDistanceAlongCurve; 
}

template <class Real_v>
inline
TemplateFieldTrack<Real_v> & TemplateFieldTrack<Real_v>::
  operator = ( const TemplateFieldTrack<Real_v>& rStVec )
{
  if (&rStVec == this) return *this;

  //try auto-vectorization
  for (int i = 0; i < 6; ++i)
  {
    fPositionMomentum[i] = rStVec.fPositionMomentum[i];
  }

  SetCurveLength( rStVec.GetCurveLength() );

  // fKineticEnergy= rStVec.fKineticEnergy;
  // fRestMass_c2= rStVec.fRestMass_c2;
  SetLabTimeOfFlight( rStVec.GetLabTimeOfFlight()  ); 
  SetProperTimeOfFlight( rStVec.GetProperTimeOfFlight()  ); 
  // SetPolarization( rStVec.GetPolarization() );
  // fMomentumDir= rStVec.fMomentumDir;

  // fCharge= rStVec.fCharge;
  return *this;
}


template <class Real_v>
TemplateFieldTrack<Real_v>::
  TemplateFieldTrack( const vecgeom::Vector3D<Real_v> & pPosition, 
                      const vecgeom::Vector3D<Real_v> & pMomentum,
                      // Real_v       restMass_c2,
                      // Real_v charge, 
                      Real_v LaboratoryTimeOfFlight,
                      Real_v curve_length          )
                      // const    ThreeVectorSIMD& vecPolarization,
                      // Real_v   curve_length )
:  fDistanceAlongCurve(curve_length),
   // fMomentumMag(pMomentum.Mag()),
   // fKineticEnergy(kineticEnergy), fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(LaboratoryTimeOfFlight), 
   fProperTimeOfFlight(0.) // ,
   // fMomentumDir(pMomentum.Unit()),
   // fCharge( charge )
{
  SetMomentum( pMomentum ); 
  SetPosition( pPosition ); 
}

// -------------------------------------------------------------------
template <class Real_v>
TemplateFieldTrack<Real_v>::TemplateFieldTrack( char )                  //  Nothing is set !!
  : // fKineticEnergy(0.), 
    // fRestMass_c2(0.), 
    fLabTimeOfFlight(0.),
    fProperTimeOfFlight(0.) // , 
    // fCharge(  DBL_MAX )
{
  vecgeom::Vector3D<Real_v> Zero(0.0, 0.0, 0.0);

  SetCurvePnt( Zero, Zero, 0.0 );
  // SetMomentum( Zero );  // Sets momentum direction as well.
  // SetPosition( Zero ); 

  // SetPolarization( Zero ); 
}

// -------------------------------------------------------------------

// Load values from array
//  
//   note that momentum direction must-be/is normalised

template <class Real_v>
void TemplateFieldTrack<Real_v>
  ::LoadFromArray(const Real_v valArrIn[ncompSVEC],
                                                    int noVarsIntegrated)
{
  int i;

  typedef vecgeom::Vector3D<Real_v> ThreeVectorSIMD;
  // Fill the variables not integrated with zero -- so it's clear !!
  // vecgeom::Vector3D<Real_v> valArr[ncompSVEC];
  Real_v valArr[ncompSVEC];
  for( i=0; i<noVarsIntegrated; i++){
     valArr[i]= valArrIn[i];
  }
  for( i=noVarsIntegrated; i<ncompSVEC; i++) {
     valArr[i]= 0.0; 
  }

#if 1  
  SetCurvePnt( ThreeVectorSIMD( valArr[0], valArr[1], valArr[2]),
               ThreeVectorSIMD( valArr[3], valArr[4], valArr[5]),
               0 ); // DistanceAlongCurve
#else  
  fPositionMomentum[0]=valArr[0];
  fPositionMomentum[1]=valArr[1];
  fPositionMomentum[2]=valArr[2];
  fPositionMomentum[3]=valArr[3];
  fPositionMomentum[4]=valArr[4];
  fPositionMomentum[5]=valArr[5];

  ThreeVectorSIMD Momentum(valArr[3],valArr[4],valArr[5]);

  // fMomentumDir= Momentum.Unit();
#endif
  
  // fKineticEnergy = momentum_square / 
  //                 (std::sqrt(momentum_square+fRestMass_c2*fRestMass_c2)
  //                  + fRestMass_c2 ); 
  // The above equation is stable for small and large momenta

  // The following components may or may not be
  //    integrated over -- integration is optional
  // fKineticEnergy= valArr[6];

  fLabTimeOfFlight=valArr[7];
  fProperTimeOfFlight=valArr[8];

  // ThreeVectorSIMD  vecPolarization= ThreeVectorSIMD(valArr[9],valArr[10],valArr[11]);
  //  SetPolarization( vecPolarization ); 

  // fMomentumDir=ThreeVectorSIMD(valArr[13],valArr[14],valArr[15]);
  // fDistanceAlongCurve= valArr[]; 
}  

template <class Real_v>
std::ostream& operator<<( std::ostream& os, const TemplateFieldTrack<Real_v>& SixVec)
{
     typedef vecgeom::Vector3D<Real_v> ThreeVectorSIMD;

     const Real_v *SixV = SixVec.fPositionMomentum;
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " "
                  << SixV[2] << " ";  // Position
     os << " P= " << SixV[3] << " " << SixV[4] << " "
                  << SixV[5] << " ";  // Momentum
     ThreeVectorSIMD  momentum(SixV[3], SixV[4], SixV[5]);
     Real_v   momentumMag= momentum.Mag();
     os << " Pmag= " << momentumMag;     
     // os << " Ekin= " << SixVec.fKineticEnergy ;
     // os << " m0= " <<   SixVec.fRestMass_c2;
     os << " Pdir= " <<  ( momentumMag > 0 ? momentum.Unit() : momentum );
     // os << " PolV= " << SixVec.GetPolarization(); 
     os << " l= " <<    SixVec.GetCurveLength();
     os << " t_lab= " << SixVec.fLabTimeOfFlight; 
     os << " t_proper= " << SixVec.fProperTimeOfFlight ; 
     os << " ) ";
     return os;
}


#endif  /* End of ifndef GUFieldTrack_HH */
