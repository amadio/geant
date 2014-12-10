//
// $Id: GUFieldTrack.cc 81175 2014-05-22 07:39:10Z gcosmo $
//
// -------------------------------------------------------------------

#include "GUFieldTrack.h"

GUFieldTrack::GUFieldTrack( const ThreeVector& pPosition, 
			          double       LaboratoryTimeOfFlight,
			    const ThreeVector& pMomentumDirection,
			          double       kineticEnergy,
			          double       restMass_c2,
                                  double       charge, 
                                  double       curve_length )
                            //     const ThreeVector& vecPolarization,
			    //     double       magnetic_dipole_moment,
                            //     double       curve_length,
                            //     double       pdgSpin )
:  fDistanceAlongCurve(curve_length),
   fKineticEnergy(kineticEnergy),
   fRestMass_c2(restMass_c2),
   fLabTimeOfFlight(LaboratoryTimeOfFlight), 
   fProperTimeOfFlight(0.),
   // fMomentumDir(pMomentumDirection),
   // fChargeState(  charge, magnetic_dipole_moment, pdgSpin )
   fCharge( charge )
   // fPDGSpin( pdgSpin )
{
  UpdateFourMomentum( kineticEnergy, pMomentumDirection ); 
      // Sets momentum direction as well.

  SetPosition( pPosition ); 
  // SetPolarization( vecPolarization ); 
}

// -------------------------------------------------------------------
   
GUFieldTrack::GUFieldTrack( char )                  //  Nothing is set !!
  : fKineticEnergy(0.), fRestMass_c2(0.), fLabTimeOfFlight(0.),
    fProperTimeOfFlight(0.), fCharge(  DBL_MAX )
{
  ThreeVector Zero(0.0, 0.0, 0.0);
  SetCurvePnt( Zero, Zero, 0.0 );
  //  SetPolarization( Zero ); 
  // fInitialMomentumMag= 0.0; // Invalid
  // fLastMomentumMag= 0.0; 
}

// -------------------------------------------------------------------

// Load values from array
//  
//   note that momentum direction must-be/is normalised

void GUFieldTrack::LoadFromArray(const double valArrIn[ncompSVEC],
                                       int noVarsIntegrated)
{
  int i;

  // Fill the variables not integrated with zero -- so it's clear !!
  double valArr[ncompSVEC];
  for( i=0; i<noVarsIntegrated; i++){
     valArr[i]= valArrIn[i];
  }
  for( i=noVarsIntegrated; i<ncompSVEC; i++) {
     valArr[i]= 0.0; 
  }

  SixVector[0]=valArr[0];
  SixVector[1]=valArr[1];
  SixVector[2]=valArr[2];
  SixVector[3]=valArr[3];
  SixVector[4]=valArr[4];
  SixVector[5]=valArr[5];

  ThreeVector Momentum(valArr[3],valArr[4],valArr[5]);

  double momentum_square= Momentum.Mag2();
  fMomentumDir= Momentum.Unit();

  fKineticEnergy = momentum_square / 
                   (std::sqrt(momentum_square+fRestMass_c2*fRestMass_c2)
                     + fRestMass_c2 ); 
  // The above equation is stable for small and large momenta

  // The following components may or may not be
  //    integrated over -- integration is optional
  // fKineticEnergy= valArr[6];

  fLabTimeOfFlight=valArr[7];
  fProperTimeOfFlight=valArr[8];
  ThreeVector  vecPolarization= ThreeVector(valArr[9],valArr[10],valArr[11]);
  //  SetPolarization( vecPolarization ); 

  // fMomentumDir=ThreeVector(valArr[13],valArr[14],valArr[15]);
  // fDistanceAlongCurve= valArr[]; 
}  

std::ostream& operator<<( std::ostream& os, const GUFieldTrack& SixVec)
{
     const double *SixV = SixVec.SixVector;
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " "
                  << SixV[2] << " ";  // Position
     os << " P= " << SixV[3] << " " << SixV[4] << " "
                  << SixV[5] << " ";  // Momentum
     os << " Pmag= "
        << ThreeVector(SixV[3], SixV[4], SixV[5]).Mag(); // mom magnitude
     os << " Ekin= " << SixVec.fKineticEnergy ;
     os << " m0= " <<   SixVec.fRestMass_c2;
     os << " Pdir= " <<  SixVec.fMomentumDir.Mag();
     os << " PolV= " << SixVec.GetPolarization(); 
     os << " l= " <<    SixVec.GetCurveLength();
     os << " t_lab= " << SixVec.fLabTimeOfFlight; 
     os << " t_proper= " << SixVec.fProperTimeOfFlight ; 
     os << " ) ";
     return os;
}
