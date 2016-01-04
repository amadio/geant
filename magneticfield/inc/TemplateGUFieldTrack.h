//
// class TemplateGUFieldTrack
//
// Class description:
//
// Data structure bringing together a magnetic track's state.
// (position, momentum direction & modulus, energy, spin, ... )
// Uses/abilities:
//  - does not maintain any relationship between its data (eg energy/momentum).
//  - for use in Runge-Kutta solver (in passing it the values right now).

// History  - derived from G4FieldTrack
// - First version: Dec 9, 2014 John Apostolakis
// -------------------------------------------------------------------

#ifndef TemplateGUFieldTrack_HH
#define TemplateGUFieldTrack_HH

#include "base/Vector3D.h"   // VecGeom/base/Vector3D.h 


// #include "G4ChargeState.hh"

template <class Backend>
class  TemplateGUFieldTrack
{

   public:  // with description

     typedef vecgeom::Vector3D<typename Backend::precision_v> ThreeVector;
     typedef typename Backend::precision_v Double_v;

     TemplateGUFieldTrack(const ThreeVector& pPosition, 
                          const ThreeVector& pMomentum,
                         // double       restMass_c2,
                          Double_v                  charge=-1,
                          Double_v  laboratoryTimeOfFlight= 0.0,
                          Double_v            curve_length= 0.0); 

     TemplateGUFieldTrack( const TemplateGUFieldTrack<Backend>&   pFieldTrack ); 
     TemplateGUFieldTrack( char );   //  Almost default constructor

     ~TemplateGUFieldTrack();
       // End of preferred Constructors / Destructor 

     inline void
     UpdateState( const ThreeVector& pPosition, 
                        Double_v     LaboratoryTimeOfFlight,
                  const ThreeVector& pMomentumDirection,
                        Double_v     momentum              ); 
        //  Update four-vectors for space/time and momentum/energy
        //    Also resets curve length.

     // void SetCharge(double charge) { fCharge= charge; }

     inline TemplateGUFieldTrack& operator = ( const TemplateGUFieldTrack<Backend> & rStVec );
       // Assignment operator

     inline ThreeVector  GetMomentum() const;
     inline ThreeVector  GetPosition() const;
     inline ThreeVector  GetMomentumDirection() const;
     inline Double_v     GetMomentumMag() const;
     inline Double_v     GetCurveLength() const;
       // Distance along curve of point.

     // inline ThreeVector  GetPolarization()   const; 
     // inline void         SetPolarization( const ThreeVector& vecPol );

     inline Double_v     GetLabTimeOfFlight() const;
     inline Double_v     GetProperTimeOfFlight() const;
     // inline double       GetKineticEnergy() const;
       // Accessors.

     inline void SetPosition(ThreeVector nPos); 
     inline void SetMomentum(ThreeVector nMom);
       // Does change mom-dir too.

     inline void SetCurvePnt(const ThreeVector& pPosition, 
                             const ThreeVector& pMomentum,  
                                   Double_v      s_curve );

     // inline void SetMomentumDir(ThreeVector nMomDir);
     // Does NOT change Momentum or Velocity Vector.

     // inline void SetRestMass(double Mass_c2) { fRestMass_c2= Mass_c2; }

       // Access
     // inline double GetCharge() const { return fCharge; } 
   
     inline void SetCurveLength(Double_v nCurve_s);
       // Distance along curve.
     inline void SetKineticEnergy(Double_v nEnergy);
       // Does not modify momentum.

     inline void SetLabTimeOfFlight(Double_v tofLab); 
     inline void SetProperTimeOfFlight(Double_v tofProper);
       //  Modifiers

   public: // without description

     static constexpr int ncompSVEC = 12;
       // Needed and should be used only for RK integration driver

     inline void DumpToArray(Double_v valArr[ncompSVEC]) const; 
            void LoadFromArray(const Double_v valArr[ncompSVEC],
                                     int      noVarsIntegrated );
     template <class Backend_>
     friend  std::ostream&
             operator<<( std::ostream& os, const TemplateGUFieldTrack<Backend_>& SixVec);

     Double_v  fPositionMomentum[6]; //initially SixVector[6]

   private:

     Double_v  fDistanceAlongCurve;  // distance along curve of point
     Double_v  fMomentumMag;
     // double  fKineticEnergy;
     // double  fRestMass_c2;
     Double_v  fLabTimeOfFlight;
     Double_v  fProperTimeOfFlight;
     // ThreeVector fPolarization;
     // ThreeVector fMomentumDir;
     // double  fInitialMomentumMag;  // At 'track' creation.
     // double  fLastMomentumMag;     // From last Update (for checking.)

     // double fCharge;


}; 

// #include "TemplateGUFieldTrack.icc"

//
// $Id: TemplateGUFieldTrack.icc 81175 2014-05-22 07:39:10Z gcosmo $
//


template <class Backend>
inline
TemplateGUFieldTrack<Backend>::
  TemplateGUFieldTrack( const TemplateGUFieldTrack<Backend>&  rStVec  )
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

  // fpChargeState= new G4ChargeState( *rStVec.fpChargeState );
  // Can share charge state only when using handles etc
  //   fpChargeState = rStVec.fpChargeState;  
}

template <class Backend>
inline
TemplateGUFieldTrack<Backend>::~TemplateGUFieldTrack()
{
  // delete fpChargeState; 
}

template <class Backend>
inline void
TemplateGUFieldTrack<Backend>::
  SetCurvePnt(const vecgeom::Vector3D<typename Backend::precision_v> &pPosition, 
              const vecgeom::Vector3D<typename Backend::precision_v> &pMomentum,  
                                      typename Backend::precision_v     s_curve )
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


template <class Backend>
inline
vecgeom::Vector3D<typename Backend::precision_v> 
TemplateGUFieldTrack<Backend>::
  GetPosition() const
{
   vecgeom::Vector3D<typename Backend::precision_v> myPosition( fPositionMomentum[0], fPositionMomentum[1], fPositionMomentum[2] );
   return myPosition;
} 

template <class Backend>
inline
vecgeom::Vector3D<typename Backend::precision_v> 
TemplateGUFieldTrack<Backend>::
  GetMomentumDirection() const 
{
   typedef vecgeom::Vector3D<typename Backend::precision_v> ThreeVector;
   double inv_mag= 1.0 / fMomentumMag;
   return inv_mag * ThreeVector( fPositionMomentum[3], fPositionMomentum[4], fPositionMomentum[5] );
} 

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  SetPosition( vecgeom::Vector3D<typename Backend::precision_v> pPosition) 
{
  //try auto-vectorization
  for (int i = 0; i < 3; ++i)
  {
    fPositionMomentum[i] = pPosition[i];
  }
} 

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  SetMomentum( vecgeom::Vector3D<typename Backend::precision_v> vMomentum) 
{

  // try auto-vectorization
  for (int i = 0; i < 3; ++i)
  {
    fPositionMomentum[i+3] = vMomentum[i];
  }

   fMomentumMag= vMomentum.Mag();
} 

template <class Backend>
inline
typename Backend::precision_v 
TemplateGUFieldTrack<Backend>::
  GetMomentumMag() const 
{
   return fMomentumMag;
} 

template <class Backend>
inline
typename Backend::precision_v  
TemplateGUFieldTrack<Backend>::
  GetCurveLength() const 
{
     return  fDistanceAlongCurve;  
}

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  SetCurveLength(typename Backend::precision_v nCurve_s)
{
     fDistanceAlongCurve= nCurve_s;  
}

// inline double TemplateGUFieldTrack<Backend>::GetKineticEnergy() const
// { return fKineticEnergy; }

// inline void TemplateGUFieldTrack<Backend>::SetKineticEnergy(double newKinEnergy)
// {  fKineticEnergy=newKinEnergy; }

// inline ThreeVector TemplateGUFieldTrack<Backend>::GetPolarization() const
// { return fPolarization; }

// inline void TemplateGUFieldTrack<Backend>::SetPolarization(const ThreeVector& vecPlz)
// { fPolarization= vecPlz; }

template <class Backend>
inline
typename Backend::precision_v 
TemplateGUFieldTrack<Backend>::
  GetLabTimeOfFlight() const
{
   return fLabTimeOfFlight;
}

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  SetLabTimeOfFlight(typename Backend::precision_v nTOF)
{
   fLabTimeOfFlight=nTOF;
}

template <class Backend>
inline
typename Backend::precision_v  
TemplateGUFieldTrack<Backend>:: 
  GetProperTimeOfFlight() const
{
   return fProperTimeOfFlight;
}

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  SetProperTimeOfFlight(typename Backend::precision_v nTOF)
{
   fProperTimeOfFlight=nTOF;
}

template <class Backend>
inline
vecgeom::Vector3D<typename Backend::precision_v> 
TemplateGUFieldTrack<Backend>::
  GetMomentum() const 
{
   return ThreeVector( fPositionMomentum[3], fPositionMomentum[4], fPositionMomentum[5] );
} 

// Dump values to array
//  
//   note that momentum direction is not saved 

template <class Backend>
inline
void TemplateGUFieldTrack<Backend>::
  DumpToArray(typename Backend::precision_v valArr[ncompSVEC] ) const
{

  //try auto-vectorization
  for (int i = 0; i < 6; ++i)
  {
    valArr[i] = fPositionMomentum[i];
  }


  ThreeVector Momentum(valArr[3],valArr[4],valArr[5]);

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

template <class Backend>
inline
TemplateGUFieldTrack<Backend> & TemplateGUFieldTrack<Backend>::
  operator = ( const TemplateGUFieldTrack<Backend>& rStVec )
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
  // (*Fpchargestate)= *(rStVec.fpChargeState);
  // fpChargeState= rStVec.fpChargeState; // Handles!!
  return *this;
}


template <class Backend>
TemplateGUFieldTrack<Backend>::
  TemplateGUFieldTrack( const vecgeom::Vector3D<typename Backend::precision_v> & pPosition, 
                        const vecgeom::Vector3D<typename Backend::precision_v> & pMomentum,
                        // double       restMass_c2,
                                                typename Backend::precision_v /*charge*/ , 
                                                typename Backend::precision_v LaboratoryTimeOfFlight,
                                                typename Backend::precision_v curve_length          )
                        // const ThreeVector& vecPolarization,
                        // double       magnetic_dipole_moment,
                        // double       curve_length,
                        // double       pdgSpin )
:  fDistanceAlongCurve(curve_length),
   // fMomentumMag(pMomentum.Mag()),
   // fKineticEnergy(kineticEnergy), fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(LaboratoryTimeOfFlight), 
   fProperTimeOfFlight(0.) // ,
   // fMomentumDir(pMomentum.Unit()),
   // fChargeState(  charge, magnetic_dipole_moment, pdgSpin )
   // fPDGSpin( pdgSpin )
   // fCharge( charge )
{
  SetMomentum( pMomentum ); 

  SetPosition( pPosition ); 
}

// -------------------------------------------------------------------
template <class Backend>
TemplateGUFieldTrack<Backend>::TemplateGUFieldTrack( char )                  //  Nothing is set !!
  : // fKineticEnergy(0.), 
    // fRestMass_c2(0.), 
    fLabTimeOfFlight(0.),
    fProperTimeOfFlight(0.) // , 
    // fCharge(  DBL_MAX )
{
  vecgeom::Vector3D<typename Backend::precision_v> Zero(0.0, 0.0, 0.0);

  SetCurvePnt( Zero, Zero, 0.0 );
  // SetMomentum( Zero );  // Sets momentum direction as well.
  // SetPosition( Zero ); 

  // SetPolarization( Zero ); 
}

// -------------------------------------------------------------------

// Load values from array
//  
//   note that momentum direction must-be/is normalised

template <class Backend>
void TemplateGUFieldTrack<Backend>
  ::LoadFromArray(const typename Backend::precision_v valArrIn[ncompSVEC],
                                                    int noVarsIntegrated)
{
  int i;

  typedef vecgeom::Vector3D<typename Backend::precision_v> ThreeVector;
  // Fill the variables not integrated with zero -- so it's clear !!
  // vecgeom::Vector3D<typename Backend::precision_v> valArr[ncompSVEC];
  typename Backend::precision_v valArr[ncompSVEC];
  for( i=0; i<noVarsIntegrated; i++){
     valArr[i]= valArrIn[i];
  }
  for( i=noVarsIntegrated; i<ncompSVEC; i++) {
     valArr[i]= 0.0; 
  }

#if 1  
  SetCurvePnt( ThreeVector( valArr[0], valArr[1], valArr[2]),
               ThreeVector( valArr[3], valArr[4], valArr[5]),
               0 ); // DistanceAlongCurve
#else  
  fPositionMomentum[0]=valArr[0];
  fPositionMomentum[1]=valArr[1];
  fPositionMomentum[2]=valArr[2];
  fPositionMomentum[3]=valArr[3];
  fPositionMomentum[4]=valArr[4];
  fPositionMomentum[5]=valArr[5];

  ThreeVector Momentum(valArr[3],valArr[4],valArr[5]);

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

  // ThreeVector  vecPolarization= ThreeVector(valArr[9],valArr[10],valArr[11]);
  //  SetPolarization( vecPolarization ); 

  // fMomentumDir=ThreeVector(valArr[13],valArr[14],valArr[15]);
  // fDistanceAlongCurve= valArr[]; 
}  

template <class Backend>
std::ostream& operator<<( std::ostream& os, const TemplateGUFieldTrack<Backend>& SixVec)
{
     typedef vecgeom::Vector3D<typename Backend::precision_v> ThreeVector;

     const typename Backend::precision_v *SixV = SixVec.fPositionMomentum;
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " "
                  << SixV[2] << " ";  // Position
     os << " P= " << SixV[3] << " " << SixV[4] << " "
                  << SixV[5] << " ";  // Momentum
     ThreeVector momentum(SixV[3], SixV[4], SixV[5]);
     typename Backend::precision_v momentumMag= momentum.Mag();
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
