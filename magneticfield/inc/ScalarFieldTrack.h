//
// class ScalarFieldTrack
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

#ifndef ScalarFieldTrack_HH
#define ScalarFieldTrack_HH

#include "base/Vector3D.h"   // VecGeom/base/Vector3D.h 
typedef vecgeom::Vector3D<double> ThreeVector;

// #include "G4ChargeState.hh"

class  ScalarFieldTrack
{
   public:  // with description

     ScalarFieldTrack( const ThreeVector& pPosition, 
                       const ThreeVector& pMomentum,
                         // double       restMass_c2,
                       double       charge,
                         // double       laboratoryTimeOfFlight= 0.0,
                       double       curve_length= 0.0); 

     ScalarFieldTrack( const ScalarFieldTrack&   pFieldTrack ); 
     ScalarFieldTrack( char );   //  Almost default constructor

     ~ScalarFieldTrack();
       // End of preferred Constructors / Destructor 

     inline void
     UpdateState( const ThreeVector& pPosition, 
                  // double       LaboratoryTimeOfFlight,
                  const ThreeVector& pMomentumDirection,
                        double       momentum); 
        //  Update four-vectors for space/time and momentum/energy
        //    Also resets curve length.

     void SetCharge(double charge) { fCharge= charge; }

     inline ScalarFieldTrack& operator = ( const ScalarFieldTrack & rStVec );
       // Assignment operator

     inline ThreeVector  GetMomentum() const;
     inline ThreeVector  GetPosition() const;
     inline ThreeVector  GetMomentumDirection() const;
     inline double       GetMomentumMag() const;
     inline double       GetCurveLength() const;
       // Distance along curve of point.

     // inline ThreeVector  GetPolarization()   const; 
     // inline void         SetPolarization( const ThreeVector& vecPol );

     // inline double       GetLabTimeOfFlight() const;
     // inline double       GetProperTimeOfFlight() const;
     // inline double       GetKineticEnergy() const;
       // Accessors.

     inline void SetPosition(const ThreeVector& nPos); 
     inline void SetMomentum(const ThreeVector& nMom);
       // Does change mom-dir too.

     inline void SetCurvePnt(const ThreeVector& pPosition, 
                             const ThreeVector& pMomentum,  
                                   double       s_curve );

     // inline void SetMomentumDir(const ThreeVector& nMomDir);
       // Does NOT change Momentum or Velocity Vector.

     // inline void SetRestMass(double Mass_c2) { fRestMass_c2= Mass_c2; }

       // Access
     inline double GetCharge() const { return fCharge; } 
   
     inline void SetCurveLength(double nCurve_s);
       // Distance along curve.
     inline void SetKineticEnergy(double nEnergy);
       // Does not modify momentum.

     // inline void SetLabTimeOfFlight(double tofLab); 
     // inline void SetProperTimeOfFlight(double tofProper);
       //  Modifiers

   public: // without description

     static constexpr int ncompSVEC = 8;
       // Needed and should be used only for RK integration driver

     inline void DumpToArray(double valArr[ncompSVEC]) const; 
     void LoadFromArray(const double valArr[ncompSVEC],
                              int noVarsIntegrated);
     friend  std::ostream&
             operator<<( std::ostream& os, const ScalarFieldTrack& SixVec);

   private:
     double  SixVector[6];
     double  fDistanceAlongCurve;  // distance along curve of point
     double  fMomentumMag;
     double  fCharge;
     // double  fKineticEnergy;
     // double  fRestMass_c2;
     // double  fLabTimeOfFlight;
     // double  fProperTimeOfFlight;
     // ThreeVector fPolarization;
     // ThreeVector fMomentumDir;
     // double  fInitialMomentumMag;  // At 'track' creation.
     // double  fLastMomentumMag;     // From last Update (for checking.)
}; 

// #include "ScalarFieldTrack.icc"

//
// $Id: ScalarFieldTrack.icc 81175 2014-05-22 07:39:10Z gcosmo $
//
// -------------------------------------------------------------------

inline
ScalarFieldTrack::ScalarFieldTrack( const ScalarFieldTrack&  rStVec  )
 : fDistanceAlongCurve( rStVec.fDistanceAlongCurve),
   fMomentumMag( rStVec.fMomentumMag ),
   // fKineticEnergy( rStVec.fKineticEnergy ),
   // fRestMass_c2( rStVec.fRestMass_c2),
   // fLabTimeOfFlight( rStVec.fLabTimeOfFlight ), 
   // fProperTimeOfFlight( rStVec.fProperTimeOfFlight ) //, 
   // fMomentumModulus( rStVec.fMomentumModulus ),
   // fPolarization( rStVec.fPolarization ), 
   // fMomentumDir( rStVec.fMomentumDir ), 
   fCharge( rStVec.fCharge )
{
  SixVector[0]= rStVec.SixVector[0];
  SixVector[1]= rStVec.SixVector[1];
  SixVector[2]= rStVec.SixVector[2];
  SixVector[3]= rStVec.SixVector[3];
  SixVector[4]= rStVec.SixVector[4];
  SixVector[5]= rStVec.SixVector[5];

  // fpChargeState= new G4ChargeState( *rStVec.fpChargeState );
  // Can share charge state only when using handles etc
  //   fpChargeState = rStVec.fpChargeState;  
}

inline
ScalarFieldTrack::~ScalarFieldTrack()
{
  // delete fpChargeState; 
}

inline void
ScalarFieldTrack::SetCurvePnt(const ThreeVector& pPosition, 
                              const ThreeVector& pMomentum,  
                                double       s_curve )
{
  SixVector[0] = pPosition.x(); 
  SixVector[1] = pPosition.y(); 
  SixVector[2] = pPosition.z(); 

  SixVector[3] = pMomentum.x(); 
  SixVector[4] = pMomentum.y(); 
  SixVector[5] = pMomentum.z(); 

  fMomentumMag= pMomentum.Mag();
  if( fMomentumMag > 0.0 )
  {
     // fMomentumDir = (1.0/fMomentumMag) * pMomentum;
  }
  fDistanceAlongCurve= s_curve;

  // return *this;
} 

inline
ThreeVector ScalarFieldTrack::GetPosition() const
{
   ThreeVector myPosition( SixVector[0], SixVector[1], SixVector[2] );
   return myPosition;
} 

inline
ThreeVector ScalarFieldTrack::GetMomentumDirection() const 
{
   double inv_mag= 1.0 / fMomentumMag;
   return inv_mag * ThreeVector( SixVector[3], SixVector[4], SixVector[5] );
} 

inline
void ScalarFieldTrack::SetPosition( const ThreeVector& pPosition) 
{
   SixVector[0] = pPosition.x(); 
   SixVector[1] = pPosition.y(); 
   SixVector[2] = pPosition.z(); 
} 

inline
void ScalarFieldTrack::SetMomentum( const ThreeVector& vMomentum) 
{
   SixVector[3] = vMomentum.x();
   SixVector[4] = vMomentum.y();
   SixVector[5] = vMomentum.z();
   fMomentumMag= vMomentum.Mag();
} 

inline
double ScalarFieldTrack::GetMomentumMag() const 
{
   return fMomentumMag;
} 

inline
double  ScalarFieldTrack::GetCurveLength() const 
{
     return  fDistanceAlongCurve;  
}

inline
void ScalarFieldTrack::SetCurveLength(double nCurve_s)
{
     fDistanceAlongCurve= nCurve_s;  
}

// inline double ScalarFieldTrack::GetKineticEnergy() const
// { return fKineticEnergy; }

// inline void ScalarFieldTrack::SetKineticEnergy(double newKinEnergy)
// {  fKineticEnergy=newKinEnergy; }

// inline ThreeVector ScalarFieldTrack::GetPolarization() const
// { return fPolarization; }

// inline void ScalarFieldTrack::SetPolarization(const ThreeVector& vecPlz)
// { fPolarization= vecPlz; }

#if 0
inline
double ScalarFieldTrack::GetLabTimeOfFlight() const
{
   return fLabTimeOfFlight;
}

inline
void ScalarFieldTrack::SetLabTimeOfFlight(double nTOF)
{
   fLabTimeOfFlight=nTOF;
}

inline
double  ScalarFieldTrack::GetProperTimeOfFlight() const
{
   return fProperTimeOfFlight;
}

inline
void ScalarFieldTrack::SetProperTimeOfFlight(double nTOF)
{
   fProperTimeOfFlight=nTOF;
}
#endif

inline
ThreeVector ScalarFieldTrack::GetMomentum() const 
{
   return ThreeVector( SixVector[3], SixVector[4], SixVector[5] );
} 

// Dump values to array
//  
//   note that momentum direction is not saved 

inline
void ScalarFieldTrack::DumpToArray(double valArr[ncompSVEC] ) const
{
  valArr[0]=SixVector[0];
  valArr[1]=SixVector[1];
  valArr[2]=SixVector[2];
  valArr[3]=SixVector[3];
  valArr[4]=SixVector[4];
  valArr[5]=SixVector[5];

  ThreeVector Momentum(valArr[3],valArr[4],valArr[5]);

  // double mass_in_Kg;
  // mass_in_Kg = fEnergy / velocity_mag_sq * (1-velocity_mag_sq/c_squared);
  // valArr[6]= mass_in_Kg;

  // The following components may or may not be integrated.
  // valArr[6]= fKineticEnergy; 

  // valArr[6]=fEnergy;  // When it is integrated over, do this ...
  // valArr[7]=fLabTimeOfFlight;
  // valArr[8]=fProperTimeOfFlight;
  // valArr[9]=fPolarization.x();
  // valArr[10]=fPolarization.y();
  // valArr[11]=fPolarization.z();
  // valArr[]=fDistanceAlongCurve; 
}

inline
ScalarFieldTrack & ScalarFieldTrack::operator = ( const ScalarFieldTrack& rStVec )
{
  if (&rStVec == this) return *this;

  SixVector[0]= rStVec.SixVector[0];
  SixVector[1]= rStVec.SixVector[1];
  SixVector[2]= rStVec.SixVector[2];
  SixVector[3]= rStVec.SixVector[3];
  SixVector[4]= rStVec.SixVector[4];
  SixVector[5]= rStVec.SixVector[5];
  SetCurveLength( rStVec.GetCurveLength() );

  // fKineticEnergy= rStVec.fKineticEnergy;
  // fRestMass_c2= rStVec.fRestMass_c2;
  // SetLabTimeOfFlight( rStVec.GetLabTimeOfFlight()  ); 
  // SetProperTimeOfFlight( rStVec.GetProperTimeOfFlight()  ); 
  // SetPolarization( rStVec.GetPolarization() );
  // fMomentumDir= rStVec.fMomentumDir;

  fCharge= rStVec.fCharge;
  // (*Fpchargestate)= *(rStVec.fpChargeState);
  // fpChargeState= rStVec.fpChargeState; // Handles!!
  return *this;
}

#if 0   
inline void 
ScalarFieldTrack::UpdateFourMomentum( double momentum_mag, 
                                  const ThreeVector& momentumDirection )
{
  // double momentum_mag  = std::sqrt(kineticEnergy*kineticEnergy
  //                       +2.0*fRestMass_c2*kineticEnergy);
  ThreeVector momentumVector= momentum_mag * momentumDirection; 

  SetMomentum( momentumVector ); 
  // SixVector[3] = momentumVector.x(); 
  // SixVector[4] = momentumVector.y(); 
  // SixVector[5] = momentumVector.z(); 

  // fMomentumDir=   momentumDirection; // Set directly to avoid inaccuracy.
  // fKineticEnergy= kineticEnergy;
}

inline void ScalarFieldTrack::UpdateState( const ThreeVector& position, 
                                       // double             laboratoryTimeOfFlight,
                                const ThreeVector& momentumDirection,
                                double             kineticEnergy
                              )
{ 
  // SetCurvePnt( position, momentumVector, s_curve=0.0);     
  SetPosition( position); 
  // fLabTimeOfFlight= laboratoryTimeOfFlight;
  fDistanceAlongCurve= 0.0;

  UpdateFourMomentum( kineticEnergy, momentumDirection); 
}
#endif

#endif  /* End of ifndef ScalarFieldTrack_HH */
